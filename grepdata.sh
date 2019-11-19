#!/bin/bash
# trajs=$(ls | cut -f2 -d"."|  tr "\n" " ")   numbers after . newline->space
for aa in "${trajs[@]}"
do
cd SNAFU-TRAJ.$aa
echo $PWD
cd ..
done

# lowest ener found, see lowerener file 901.48626117767878
if [[ -e ener_*.dat ]];then
rm ener_*.dat en_*.dat
fi
for ((i=1;i<=702;i++))do
if [[ -e fomo_$i.inp.out ]];then
echo "fomo_$i.inp.out"
for ((st=1;st<=8;st++))do

grep "Doublet state  $st energy" fomo_$i.inp.out |  awk '{print 27.21139*($5+901.48626117767878)}' >> en_$st.dat
done

fi

done
for ((st=1;st<=8;st++))do
paste dists.dat en_$st.dat >ener_$st.dat 

done

rm en_*.dat
