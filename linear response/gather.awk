/^ total charge/{
chargeline=NR
}
NR==chargeline+2 && /^# of ion/{
flag=1
}
NR==chargeline+4 && flag==1 {
printf "%s",$4
exit
}
