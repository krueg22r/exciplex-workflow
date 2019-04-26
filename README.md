Finding the right DFT functional is difficult. In many cases, particularly 
for problems involving excited states, range-separated functionals are the way 
to go. Performance can be improved using a non-empirical tuning process to help 
the functional come closer to satisfying the DFT version of Koopman's theorem by 
adjusting the range-separation parameter. This parameter mu determines the speed 
of switching between local and exact exchange. 
This process requires two calculations each for several values of mu. A fit of 
the errors from each yields the optimal mu value for the geometry. 

This all adds up to a lot of calculations. However, this set of scripts allows easy 
calculation setup. The steps are: 

1) Collect your molecule geometries in xyz files in a directory.
2) Run rangeSep.py to set up tuning runs. 
3) Use the shell scripts generated to submit tuning runs one molecule
at a time. 
4) Once those calculations are done, use s1opt.py to set up excited-state
optimization runs. 
5) Submit the batch_submit shell script to add opt jobs to the queue. 
6) If any opt runs out of time, go into the molecule's S1opt directory 
and run optRestart.py. Submit the new .sh file (Restart is appended to the 
original job name). 
7) Extract some useful information from the calculations using dimerData.py. 
8) Enjoy your data! 

This library assumes a Maui/Torque scheduler and uses the Orca quantum 
chemistry code, but may easily be adapted for other applications. 
