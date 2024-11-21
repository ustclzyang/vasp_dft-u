NR==FNR && NR==3{
a=$4
}
NR!=FNR && FNR==3{
b=$4
printf "%s %s %s %d",a,b,(b-a),(b>a)
}
