#!/bin/bash
#An example for MPI job.
#SBATCH -J vasp
#SBATCH -o job-%j.log
#SBATCH -e job-%j.err
#SBATCH -n 96
#SBATCH --ntasks-per-node=48
#SBATCH -p normal

echo Time is `date`
echo Directory is $PWD
echo This job runs on the following nodes:
echo $SLURM_JOB_NODELIST
echo This job has allocated $SLURM_JOB_CPUS_PER_NODE cpu cores.
module add vasp/5.4.4-intel2019
#
# The real work starts here
#

cat > INCAR.DFT <<!
Global Parameters
ISTART =  1            (Read existing wavefunction, if there)
ISPIN  =  1            (Non-Spin polarised DFT)
# ICHARG =  11         (Non-self-consistent: GGA/LDA band structures)
LREAL  = .FALSE.       (Projection operators: automatic)
ENCUT  =  500        (Cut-off energy for plane wave basis set, in eV)
# PREC   =  Accurate   (Precision level: Normal or Accurate, set Accurate when perform structure lattice relaxation calculation)
LWAVE  = .TRUE.        (Write WAVECAR or not)
LCHARG = .TRUE.        (Write CHGCAR or not)
ADDGRID= .TRUE.        (Increase grid, helps GGA convergence)
LASPH  = .TRUE.        (Give more accurate total energies and band structure calculations)
PREC   = Accurate      (Accurate strictly avoids any aliasing or wrap around errors)
# LVTOT  = .TRUE.      (Write total electrostatic potential into LOCPOT or not)
# LVHAR  = .TRUE.      (Write ionic + Hartree electrostatic potential into LOCPOT or not)
# NELECT =             (No. of electrons: charged cells, be careful)
# LPLANE = .TRUE.      (Real space distribution, supercells)
# NWRITE = 2           (Medium-level output)
# KPAR   = 2           (Divides k-grid into separate groups)
# NGXF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGYF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGZF    = 300        (FFT grid mesh density for nice charge/potential plots)
NCORE = 8
LMAXMIX      = 4
 
Static Calculation
ISMEAR =  0            (gaussian smearing method)
SIGMA  =  0.05         (please check the width of the smearing)
LORBIT =  11           (PAW radii for projected DOS)
NEDOS  =  2001         (DOSCAR points)
NELM   =  60           (Max electronic SCF steps)
EDIFF  =  1E-08        (SCF energy convergence, in eV)
 
DFT-D3 Correction
IVDW   =  11           (DFT-D3 method of method with no damping)
!

v=7.00
band_diff=0.01
tmp_diff=2.0
maxn=10
a=6.00
b=7.00
flag=1 # 若b是一个确定的上界，这里需要设置为1
n=0

until [ $(echo "scale=4; ($band_diff>$tmp_diff) && (-$band_diff<$tmp_diff)" | bc) -ne 0 -o $n -ge $maxn ] # |band_diff|>|tmp_diff| 或者到达 nmax
do

cp KPOINTS.scf KPOINTS
cp INCAR.DFT INCAR
rm WAVECAR CHGCAR

cat >> INCAR <<!
LDAU         = .TRUE.
LDAUTYPE     =  2
LDAUL        =  2 -1 -1 
LDAUU        =  $v 0.00 0.00
!

# scf
mpirun vasp_std | tee runlog 

cp OUTCAR  OUTCAR.U=$v
rm WAVECAR
cp KPATH.in KPOINTS
cp INCAR.DFT INCAR
cat >> INCAR <<!
ICHARG       = 11

LDAU         = .TRUE.
LDAUTYPE     =  2
LDAUL        =  2 -1 -1 
LDAUU        =  $v 0.00 0.00
!
# nonsc
mpirun vasp_std | tee runlog 
vaspkit -task 211
cp OUTCAR  OUTCAR.U=$v.ICHARG=11
rm BAND.dat

read var1 var2 var3 var4 < <(awk -f bandgap.awk BAND_GAP_compare BAND_GAP)
echo "$v $var1 $var2 $var3 $var4" >> mylog
tmp_diff=$var3
if [ $var4 -eq 1 ]
then
b=$v
v=$(echo "scale=4; ($a+$b)/2" | bc)
flag=1
elif [ $flag -eq 0 ]
then
a=$v
b=$(echo "scale=4; $b+4.0" | bc)
v=$b
else
a=$v
v=$(echo "scale=4; ($a+$b)/2" | bc)
fi
let n++
done
