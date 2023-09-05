############################
# Conduct Gillespie simulations of a stochastic version of
# the competitive Lotka-Volterra equation for finite pops
# To illustrate noise-induced selection
# This is part of the manuscript by Bhat and Guttal 2023
# Written by Shikhara Bhat
# IISER Pune, India
# Date: Sat Jun 24 09:17:43 2023
###########################


library(GillespieSSA2) #for stochastic simulations
library(dplyr) #for manipulating dataframes
library(ggplot2) #for plotting

#change parameters as appropriate here 
model_params <- c(
  
  
  #Intrinsic birth and death rates  
  b1 = 1, #intrinsic per capita birth rate of type 1 inds
  b2 = 1, #intrinsic per capita birth rate of type 2 inds
  
  d1 = 1, #per capita death rate due to intraspecific competition of type 1 inds
  d2 = 1, #per capita death rate due to intraspecific competition in type 2 inds

  #natural selection effect
  epsilon = 0.0005, #competitive advantage to species 2
  
  #mutation rate
  mu = 0.001, #mutation rate (assumed symmetric) from one type to the other
  
  #Carrying capacity
  K = 500 #System-size parameter
)


#Define total birth and death rates
reactions <- list(
  
  #                     rate                                   effect

  
  #births
  #       intrinsic            interactions         mutation
  reaction("b1*N1   -     (1+epsilon)*N1*N2/K   +    mu*N2",   c(N1 = +1)),
  reaction("b2*N2                               +    mu*N1",   c(N2 = +1)),
  
  #deaths
  #         intrinsic         interactions        
  reaction("d1*N1*N1/K                      ",                 c(N1 = -1)),
  reaction("d2*N2*N2/K   +     N1*N2/K      ",                 c(N2 = -1))
)

##########################
#In terms of predictions, we have, for epsilon > 0:
#classical selection (increased fitness): species 2 favored over species 1
#noise-induced selection (reduced turnover): species 1 favored over species 2
########################

#Initial conditions
#Start with equal no. of inds of each type
initial_state <- c(
  N1 = round(as.integer(model_params['K'])/2),
  N2 = round(as.integer(model_params['K'])/2)
)


final_time <- 1e5  #time at which the Gillespie simulation stops
runs <- 100         #no. of independent realizations to simulate


#Run the SSA

set.seed(2878) #for reproducibility

start_time = Sys.time() #for timing how long running the code takes


#Time complexity: This should scale linearly with time but exponentially with K :(
sim_run <- ssa(
  initial_state = initial_state,
  reactions = reactions,
  params = model_params,
  final_time = final_time,
  method = ssa_exact(),
  census_interval = 50, #how often to store data. 0 means store everything.
  sim_name = 'Example model from Bhat and Guttal 2023'
)

#record realization number (for plotting timeseries)
sim_run_num_col <- rep(1,length(sim_run[['time']]))

sim_data <- as.data.frame(cbind(sim_run_num_col,sim_run[['time']], sim_run[['state']]))
rm(sim_run)#to save memory

#loop to run multiple realizations
for (i in 1:runs-1){
  sim_run <- ssa(
    initial_state = initial_state,
    reactions = reactions,
    params = model_params,
    final_time = final_time,
    method = ssa_exact(),
    census_interval = 50, #how often to store data. 0 means store everything.
    sim_name = 'Example model from Bhat and Guttal 2023'
  )
  
  #record realization number (for plotting timeseries)
  sim_run_num_col <- rep(i+1,length(sim_run[['time']])) 
  
  sim_data <- rbind(sim_data,as.data.frame(cbind(sim_run_num_col,sim_run[['time']], sim_run[['state']])))
  
  rm(sim_run)#to save memory
}

#Report how much time it took to run the simulations
running_time = Sys.time() - start_time
running_time


#Collate to dataframe

colnames(sim_data) <- c('run','time','N1','N2')

#move from (N1,N2) space to (p, N1+N2) space
sim_data %>% 
  mutate(Ntot = N1+N2) %>% #total population size
  mutate(p = N1/(Ntot)) %>% #proportion of type 1 individuals
  select(run,time,p,Ntot) -> transformed_data

#save to file
write.csv(transformed_data,paste('./data/K_',as.integer(model_params['K']),'_epsilon_',as.numeric(model_params['epsilon']),'_time_',final_time,'_num_runs_',runs,'.csv',sep=''),row.names=FALSE)

