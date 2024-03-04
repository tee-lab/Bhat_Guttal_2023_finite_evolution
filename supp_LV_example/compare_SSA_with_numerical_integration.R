############################
# Plot time series and steady state densities of the simulation output
# Compare IbM results with predictions obtaiend by numerically solving the
# associated Fokker-Planck equation under a separation of timescales assumption
# To illustrate noise-induced selection
# input is a CSV file obtained from noise_induced_selection_model_SSA.R
# This is part of the manuscript by Bhat and Guttal 2023
# Written by Shikhara Bhat
# IISER Pune, India
# Date: Mon Jun 26 13:37:51 2023
###########################


library(ggplot2) #for plotting
library(dplyr)   #for handling data


#Import the data
#change parameters as appropriate here
K = 500
final_time = 1e5
eps = 0.0005
mu = 0.001
runs = 100

data <- read.csv(paste('./data/K_',as.integer(K),'_epsilon_',eps,'_time_',final_time,'_num_runs_',runs,'.csv',sep=''))

#convert realization index to factor so ggplot interprets it correctly
data$run <- as.factor(data$run)

data %>%
  filter(Ntot > 0) %>%  #condition on non-extinction
  filter(time > 200) %>%  #remove initial transient
  select(time,p,Ntot,run) -> data


#get histogram of the IbM data (for plotting)
num_bins = 1000 #number of bins for the histogram
p_bins <- seq(0,1,length.out=num_bins)

hist <- hist(data$p, breaks = p_bins, plot = FALSE)
IbM_data <- data.frame(p_IbM = hist$mids, m_IbM = hist$density/sum(hist$density))

#####################################
#Numerical integration


f <- function(p){ #drift
  p*((1-p)**2)*((2-eps*(K-1))/K)+mu*(1-2*p)*(1-(1/K))
}

g <- function(p){ #square of diffusion
  (p*(1-p)*(2+(2+eps)*(-(p**2) + 2*p - 1))+mu*(1-3*p*(1-p)))/K
}

dU <- function(p){
  f(p)/g(p)
}

U <- function(p){ #potential
  integrate(Vectorize(dU),lower=0,upper=p)$value
}

m <- function(p){ #speed density
  exp(U(p))/g(p)
}


#compute numeric values
m_numeric <- c(0)

#we are integrating over the same intervals used by the histogram for binning
for (i in 2:length(p_bins)){
  
  m_numeric <- c(m_numeric,integrate(Vectorize(m), lower = p_bins[i-1], upper = p_bins[i])$value)
}


#normalize
m_numeric <- m_numeric/sum(m_numeric)

#numeric_data <- data.frame(p_numeric = IbM_data$p_IbM, m_numeric = m_numeric)
numeric_data <- data.frame(p_numeric = c(0,IbM_data$p_IbM), m_numeric = m_numeric)



###########################
#plotting
#Plot of probability density function



#only plot some points of IbM data, for clarity
 IbM_data %>%
   filter((row_number() %% 8 == 1) | p_IbM > 0.99) -> IbM_data


#plotting

dens <- ggplot(data = IbM_data) + geom_point(aes(x=p_IbM,y=m_IbM, color = 'Exact IbM (simulation)'),size = 4, alpha=0.7)
dens <- dens + geom_line(data = numeric_data, aes(x=p_numeric,y=m_numeric, color = 'Solution to approximate FPE (Theory)'),linewidth=1, linetype = 'dashed')


#aesthetics
dens <- dens + theme_light() + xlab('Proportion of Type 1 individuals') + ylab("Probability Density")
dens <- dens + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dens <- dens + theme(axis.text = element_text(face = 'bold', color = 'black',size = 30))
dens <- dens + theme(axis.title = element_text(face = 'bold',color = 'black',size = 35))
dens <- dens + scale_color_manual(values =  c('#0000ff','#ff0000'))
#dens <- dens + theme(legend.position = 'none') #remove legend
dens <- dens + theme(aspect.ratio = 1) #set aspect ratio
dens <- dens + theme(axis.ticks.length = unit(.15, "cm")) + theme(axis.ticks = element_line(color='black',linewidth=0.5))
dens <- dens + scale_y_continuous(expand = c(0, 0)) #this removes the gap between plot and x axis
dens <- dens + scale_x_continuous(limits = c(0,1),breaks=seq(0,1,by=0.1))




#Note down the infinite population limit solution
#This solution is found symbolically in the Jupyter Notebook analytical_solution_in_deterministic_limit.ipynb


#det_soln <- 0.18912081520889 #This is for eps = 0.005, mu = 0.001
det_soln <- 0.466823165069081 #This is for eps = 0.0005, mu = 0.001


#add a vertical line at the infinite population limit solution
dens <- dens + geom_vline(xintercept = det_soln,linetype='solid',linewidth=1.1,alpha=1,color='black')


#set the scale
dens <- dens + ylim(0,0.04) #for K = 500, this is the right scale 
#dens <- dens + ylim(0,0.01) #for K = 5000, this is the right scale

#save to file
ggsave(paste('./plots/K_',as.integer(K),'_epsilon_',eps,'_time_',final_time,'.svg',sep=''),dens,width=12,height=12,dpi=600)


