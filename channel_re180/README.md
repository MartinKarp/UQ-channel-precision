This computes a channel for Re_tau=180. We restart the simulation from a checkpoint and compute the time series for 9 different elements together with the statistics.

PLease note that the time series is a hack.

The time series is output directly in the log file so the time seres for a certain point can be extracted with awk e.g.

awk '/global-element-id/ {getline;if(1 == $1){print}}' logfile.out

which outputs the interesting quantities for element 1
