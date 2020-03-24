name=$1
top=$2
traj=$3
mkdir $name
cpptraj $top<<EOF
trajin $traj
strip !(@C1') 
trajout $name/C.mdcrd 
EOF
