this restart the simulation from a checkpoint and computes the time series for 9 different elements
we also compute the statstics etc., but the chacky part is the time series.
the time series is output directly in the log file so the time seres for a certain point can be extracted with awk e.g.
awk '/global-element-id/ {getline;if(1 == $1){print}}' logfile.out
which outputs the interesting qunantites for element 1
