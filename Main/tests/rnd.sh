# shell file for running tests on randomly-generated (CQKnP) instances
# nruns mx_size mn_size chgprc sr

set sr = 1
# sorting procedure, default QS, change into sr = 0 for BS

foreach mx_size ( 100000 10000 1000 100 )
foreach chgprc ( 1 0.1 0.05 0.01 ) 

switch ( $mx_size )
 case 100000:
  @ nruns = 100
  breaksw
 case 100: 
  @ nruns = 10000
  breaksw
 default:
  @ nruns = 1000
endsw

@ mn_size= $mx_size / 100

./CQKnPSolve $nruns $mx_size $mn_size $chgprc $sr

end
end

