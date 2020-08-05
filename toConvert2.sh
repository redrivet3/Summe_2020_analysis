#!/bin/csh -f
#source /group/clas12/packages/setup.csh
#module avail
#module load coatjava/6.5.3
#module load clas12root/1.4
#module load root/6.14.04
#module load groovy/2.5.6
source /group/clas12/packages/setup.csh
module avail
module load clas12/pro
module list
clas12root -q simpleAnaLC.C+ --in=/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v0/dst/recon/005038/rec_clas_005038.evio.00485-00489.hipo



