# !/bin/bash

DATEPLUS=$(date +"%Y%m%d%H%M%S.%N")

function slaptime(){
	FILENAME=test/$1.$DATEPLUS".csv"
	echo "" > $FILENAME

	echo n_loop,n_files,n_threads,t_rdfiles,s_training,e_training,s_prediction,e_prediction,rmse,t_total >> $FILENAME

	for t in 1 2 4 8; do # training period in years
		for i in $(seq 1 $2); do # threads number
			for j in $(seq 1 $3); do # how many times
			echo "countdown - test" $i - $j
			sleep 5

			echo $j,$(bin/generic_app $i $t $(ls support/nc_data/-*.nc)) >> $FILENAME

			done
		done
	done
}

$1 $2 $3 $4
# echo 0:$0 1:$1 2:$2 3:$3 4:$4 5:$5