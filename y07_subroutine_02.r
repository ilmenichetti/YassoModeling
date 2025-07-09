# In this file a light R-function of Yasso07 is created

# Copyright (C) <2017>  <Finnish Meteorological Institute>
#   
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


# Based on Yasso07 description Tuomi & Liski 17.3.2008
# Created by Taru Palosuo in December 2011


#  Instructions  IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

# 1) first run the source code for the function with "source(....r)"
# 2) then you can use the yasso07-function by just calling it yasso07.light(..)
# 3) Input needed for the function is the same as for the matrix version:
#        1. MeanTemperature - mean annual temperatures [C]
#        2. TemperatureAmplitude - temperature amplitudes [C] 
#        3. Precipitation - annual precipiations [mm]
#        4. InitialCPool - initial C pools of model compartments, length 5, [any unit]
#        5. LitterInput - mean litter input, 5 columns AWENH, [any unit]
#        6. WoodySize - size of woody litter (for non-woody litter this is 0)
#        7. Yasso07Parameters - these in the format applied in the fortran version, length 44
#        8. Time - simulation time 

# NOTE that this function eats only one type of material at the time. So, non-woody and different woody litter
# materials needs to be calculated separately.

# The output of the function is the matrix AWENH compartments at the given time since the simulation start


# Basics  BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

# additional R libraries (as needed)
library(Matrix)


#testing Git changes

# Function definition   FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

# The input for this function is provided as vectors, row = year
N=50
MeanTemperature = rep(20, N)
TemperatureAmplitude = rep(5, N)
Precipitation = rep(0.6, N)
InitialCPool = c(0,0,0,0,0)
LitterInput = rep(3, N)
WoodySize = 4.5
Yasso07Parameters = Yasso07Parameters_load$value
SimulationTime = seq(1:N)

# this function modifies the initial yasso07.light to run iteratively for each time step
yasso07.iterative = function(MeanTemperature, TemperatureAmplitude, Precipitation, InitialCPool, LitterInput, WoodySize, Yasso07Parameters, SimulationTime) {
      
      # Ensure input sequences are of the same length
      if (length(MeanTemperature) != length(SimulationTime) ||
          length(TemperatureAmplitude) != length(SimulationTime) ||
          length(Precipitation) != length(SimulationTime) ||
          length(LitterInput) != length(SimulationTime)) {
            stop("All input sequences must be of the same length.")
      }
      
      PA = Yasso07Parameters
      WS = WoodySize
      num_compartments = length(InitialCPool)
      num_time_points = length(SimulationTime)
      
      # Initialize a 3-dimensional array to store the results
      LC_array = array(0, dim = c(num_compartments, num_time_points))
      
      alfa = c(-PA[1], -PA[2], -PA[3], -PA[4], -PA[35])  # Vector of decomposition rates
      # Creating the matrix A_p (here called p)
      row1 = c(-1, PA[5], PA[6], PA[7], 0)
      row2 = c(PA[8], -1, PA[9], PA[10], 0)
      row3 = c(PA[11], PA[12], -1, PA[13], 0)
      row4 = c(PA[14], PA[15], PA[16], -1, 0)
      row5 = c(PA[36], PA[36], PA[36], PA[36], -1)
      
      p = matrix(c(row1, row2, row3, row4, row5), 5, 5, byrow = TRUE)
      
      # Temperature and moisture dependence parameters
      beta1 = PA[17]
      beta2 = PA[18]
      gamma = PA[26]
      
      # Woody litter size dependence parameters
      delta1 = PA[39]
      delta2 = PA[40]
      r = PA[41]
      
      # Initialize the initial carbon pool
      LC = matrix(0, nrow = num_compartments, ncol = num_time_points)
      LC[,1]= InitialCPool
      
      for (i in 2:length(SimulationTime)) {
         
         #reding the transient variables for time step i
            MT = MeanTemperature[i]
            TA = TemperatureAmplitude[i]
            PR = Precipitation[i] / 1000  # conversion from mm to meters
            Previous_state =  LC[,i-1] #the state to be updated
            LI = LitterInput[i]
            TI = SimulationTime[i]
            
            
            T1 = MT + 4 * TA / pi * (1 / sqrt(2) - 1)          # Eq. 2.4 in model description
            T2 = MT - 4 * TA / (sqrt(2) * pi)                  # Eq. 2.5 in model description
            T3 = MT + 4 * TA / pi * (1 - 1 / sqrt(2))          # Eq. 2.6 in model description
            T4 = MT + 4 * TA / (sqrt(2) * pi)                  # Eq. 2.7 in model description
            
            # k following Eq. 3 in Tuomi et al. 2009. Eco.Mod. 220: 3362-3371
            k = alfa * mean(exp(beta1 * c(T1, T2, T3, T4) + beta2 * (c(T1, T2, T3, T4)^2)) * (1 - exp(gamma * PR)))
            
            # The effect of wl size as in Eq. 3.1 in model description
            k = c(k[1:4] * (1 + delta1 * WS + delta2 * (WS^2))^(r), k[5])
            
            A = p %*% diag(k)  # Matrix multiplication in R: %*%
            
            #this runs for one time step only
            TI=1
            
            # Analytical solution as in Eq. 1.3 in model description
            LC_step = solve(A) %*% (expm(A * TI) %*% (A %*% Previous_state + LI) - LI)
            
            # Store the result in the 3-dimensional array
            LC[, i] = as.vector(LC_step)
      }
      
      return(LC)
}


