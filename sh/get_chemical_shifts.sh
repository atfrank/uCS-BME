#!/bin/bashx
if [[ $# -ne 5 ]]
then
	echo "usage: $0 <id>  <#models> <path-to-ct> <ss2cs-path> <output>"
else
	id=$1
	models=`seq 1 ${2}`	
	ct_data=$3	
	SS2CS=$4
	output=$5
	
	# initialize new file
	prefix=${ct_data}/${id}
	echo "model,resid,resname,nucleus,simcs,id" > ${output}
	
	# loop over models and compute chemical shifts
	for model in ${models}
	do
		if [[ -f ${prefix}_${model}.ct ]]
		then
			# get chemical shifts
			python ${SS2CS}/ss2cs.py ${prefix}_${model}.ct ${id} ${prefix}.csv ${model} ${SS2CS} &> /dev/null
			cat ${prefix}.csv >> ${output}
		fi
	done
fi