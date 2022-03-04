#!/bin/bashx
if [[ $# -ne 4 ]]
then
	echo "usage: $0 <id>  <#models> <path-to-ct> <ss2cs-path> "
else
	id=$1
	models=`seq 1 ${2}`	
	ct_data=$3	
	SS2CS=$4
	
	# initialize new file
	prefix=${ct_data}/${id}
	echo "model,resid,resname,nucleus,simcs,id" > ${prefix}_formatted.csv	

	# loop over models and compute chemical shifts
	for model in ${models}
	do
		if [[ -f ${prefix}_${model}.ct ]]
		then
			# get chemical shifts
			python ${SS2CS}/ss2cs.py ${prefix}_${model}.ct ${id} ${prefix}.csv ${model} ${SS2CS} &> /dev/null
			cat ${prefix}.csv >> ${prefix}_formatted.csv		
		fi
	done
	head ${prefix}_formatted.csv	
fi