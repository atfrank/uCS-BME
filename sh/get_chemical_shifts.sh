#!/bin/bashx
if [[ $# -ne 4 ]]
then
	echo "usage: $0 <id>  <#models> <path-to-ct> <ss2cs-path> "
else
	id=$1
	models=`seq 1 ${2}`	
	ct_data=$3	
	SS2CS=$4
	
	# clean up
	rm -f ${prefix}_formatted.csv		
	
	for model in ${models}
	do
		prefix=${ct_data}/${id}
		if [[ -f ${prefix}_${model}.ct ]]
		then
			# get chemical shifts
			python ${SS2CS}/ss2cs.py ${prefix}_${model}.ct ${id} ${prefix}.csv  ${SS2CS} &> /dev/null
			awk -v model=$model  '{print model, $0}' ${prefix}.csv | tee -a ${prefix}_formatted.csv		
		fi
	done
fi