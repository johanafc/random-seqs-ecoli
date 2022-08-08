rm -v ids*

while read -r line;
do
  echo $line | awk '{print ">"$1"\n"$8}' >test
  echo "Running IUPred -short"
  python3 ~/software/iupred/iupred2a.py test short>test2
  grep -v "#" test2| awk '{ sum += $3 } END { if (NR > 0) print sum / NR }'>>ids-short-avg
  grep -v "#" test2| awk 'BEGIN {a=0} {if ($3 >= 0.5) a++} END { if (NR > 0) print a / NR }'>>ids-short-frac

  echo "Running IUPred -long"
  python3 ~/software/iupred/iupred2a.py test long>test3
  grep -v "#" test3| awk '{ sum += $3 } END { if (NR > 0) print sum / NR }'>>ids-long-avg
  grep -v "#" test3| awk 'BEGIN {a=0} {if ($3 >= 0.5) a++} END { if (NR > 0) print a / NR }'>>ids-long-frac
done<$1

rm -v test*
