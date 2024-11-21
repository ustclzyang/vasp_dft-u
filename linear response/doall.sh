#!/bin/bash
#An example for MPI job.
#SBATCH -J vasp
#SBATCH -o job-%j.log
#SBATCH -e job-%j.err
#SBATCH -n 96
#SBATCH --ntasks-per-node=48
#SBATCH -p short

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
ISPIN  =  2            (Non-Spin polarised DFT)
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

cp INCAR.DFT INCAR
rm WAVECAR CHGCAR

mpirun vasp_std | tee runlog 

cp OUTCAR  OUTCAR.0
cp OSZICAR OSZICAR.0
cp WAVECAR WAVECAR.0
cp CHGCAR  CHGCAR.0


for v in +0.05 -0.05 +0.10 -0.10 +0.15 -0.15 +0.20 -0.20
do

cp INCAR.DFT INCAR
cat >> INCAR <<!
ICHARG       = 11

LDAU         = .TRUE.
LDAUTYPE     =  3
LDAUL        =  2 -1 -1 -1
LDAUU        =  $v 0.00 0.00 0.00
LDAUJ        =  $v 0.00 0.00 0.00
LDAUPRINT    =  2
!

cp WAVECAR.0 WAVECAR
cp CHGCAR.0  CHGCAR

mpirun vasp_std | tee runlog 

cp OSZICAR OSZICAR.V=$v.ICHARG=11
cp OUTCAR  OUTCAR.V=$v.ICHARG=11

cp INCAR.DFT INCAR
cat >> INCAR <<!
LDAU         = .TRUE.
LDAUTYPE     =  3
LDAUL        =  2 -1 -1 -1
LDAUU        =  $v 0.00 0.00 0.00
LDAUJ        =  $v 0.00 0.00 0.00
LDAUPRINT    =  2
!

mpirun vasp_std | tee runlog 

cp OSZICAR OSZICAR.V=$v
cp OUTCAR  OUTCAR.V=$v

done
