TOTAL_TASKS=1
PREFIX=sample

# cd ..
EXP_LIST=./outputs/${PREFIX}*
IFS=$'\n' sorted=($(sort -r <<<"${EXP_LIST[*]}"))

if [ $# == 3 ]; then
    NODE_ALL=$1
    NODE_THIS=$2
    START_IDX=$3
fi
if [ $# == 2 ]; then
    NODE_ALL=$1
    NODE_THIS=$2
    START_IDX=0
fi

#NODE_ALL=1
#NODE_THIS=2
#START_IDX=0

#for ((i=$START_IDX;i<$TOTAL_TASKS;i++)); do
for ((i=START_IDX;i<TOTAL_TASKS;i++)); do
    echo $NODE_ALL, $i
    NODE_TARGET=$(($i % $NODE_ALL))
    if [ $NODE_TARGET == $NODE_THIS ]; then
        exprdir="${sorted[i]}"
        echo $exprdir
        exprname=`basename ${exprdir}`
        echo $exprname
        if [ ! -f "${exprdir}/results.pt" ]; then
            python evaluation/evaluate_chem.py --exp_name ${exprname}
        fi
    fi
done
