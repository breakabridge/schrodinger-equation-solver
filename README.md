# schrodinger-equation-solver
A 2D C++ solver of Schrödinger's equation for a wavefunction ψ in a scalar potential V. The code uses the FTCS method to calculate ψ at each timestep. There is Python code to animate |ψ|². The idea comes from the 2011 paper "Numerical Resolution of the Schrödinger Equation" by Jørgensen et al. Includes example videos I created.

How to use:
  -- Compile and run "schrodinger.cpp". You may wish to change the simulation parameters as described in the code comments. This will output a data file,     
      "schrodinger_sim.dat".
      
  -- Run "schrodinger_plotter.py". This reads the data file and creates the video. NOTE: THE PROGRAM Ffmpeg MUST BE INSTALLED FOR THIS TO WORK. If you know some 
      matplotlib you may want to change the code so that it outputs to a .gif, etc.
