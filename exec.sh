inc=20
for i in {2..23..${inc}}; do
    num1=$(( ${i} + 1 ))
    num2=$(( ${i} + ${inc} ))

    cut preproc_tpm.gct -f 1,2,${num1}-${num2} > ./intermed/${num1}-${num2}.tsv
done
