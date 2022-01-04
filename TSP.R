
# TSP
rm(list = ls(all=T))
# load data
y <- as.matrix(eurodist) 
# defining functions
# get intital solution
get_initial_x <- function(ncity){ 
  tour <- sample(1:ncity)
  # ensure tour returns to first city
  tour <- c(tour, tour[1]) 
  return(tour)
}
cur_tour <- get_initial_x(nrow(y)) 
cur_tour

# evaluation function
evaluate_x <- function(dists, tour){ 
  sum(dists[cbind(tour[-length(tour)], tour[-1])])
}
evaluate_x(dists = y, tour = cur_tour)

# perturbation functions
# reverse operator
perturb_x <- function(tour){
  # select two cities at random
  ch <- sort(sample(2:(length(tour)-1),2,replace=FALSE)) 
  # reconnect other way
  tour[ch[1]:ch[2]] <- tour[ch[2]:ch[1]]
  return(tour) 
}

#############################
### Simulation
# Case 1.1 (base case)
# set start temperature and geometric cooling factor 
set.seed(100) # for repeatability
start_temp <- 50000
temp_factor <- 0.98
# get an initial solution
cur_x <- get_initial_x(ncol(y))
# evaluate the solution
cur_fx <- evaluate_x(dists = y, tour = cur_x)
# initialize results data frames
all_fx <- c()
all_x <- data.frame()
all_temp <- c()

# for a fixed number of iterations
set.seed(2020)
for(i in 1:10000){
  # generate a candidate solution
  prop_x <- perturb_x(cur_x)
  # evaluate the candidate solution
  prop_fx <- evaluate_x(dists = y, tour = prop_x)
  # calculate the probability of accepting the candidate
  anneal_temp <- start_temp * temp_factor ^ i 
  accept_prob <- exp(-(prop_fx - cur_fx) / anneal_temp)
  # accept or reject the candidate
  if(prop_fx < cur_fx){ 
    cur_x <- prop_x 
    cur_fx <- prop_fx
  }
  else{ if(runif(1) < accept_prob){
    cur_x <- prop_x
    cur_fx <- prop_fx
    }
  }
  # store all results
  all_fx <- c(all_fx, cur_fx) 
  all_x <- rbind(all_x,cur_x)
  all_temp <- c(all_temp, anneal_temp)
}

# plot
pdf(file = "TSPCase11.pdf", width = 7, height = 7, compress = F)
plot(all_fx, type="l", lwd=2, col="darkred", 
     xlab="Iteration", ylab="Distance") # converges quickly 
dev.off()
plot(all_temp, type="l")
# results 
best_solution <- which.min(all_fx) 
best_tour <- all_x[best_solution,] 
optimal_solution <- all_fx[best_solution] 
optimal_solution #12 919
best_tour
colnames(y)[as.numeric(best_tour)]
###############################
## Case 1.2
# explore feature space more
# increase cooling factor 
# set start temperature and geometric cooling factor 
set.seed(100) # for repeatability
start_temp <- 50000
temp_factor <- 0.999
# get an initial solution
cur_x <- get_initial_x(ncol(y))
# evaluate the solution
cur_fx <- evaluate_x(dists = y, tour = cur_x)
# initialize results data frames
all_fx <- c()
all_x <- data.frame()
all_temp <- c()

# for a fixed number of iterations
set.seed(2020)
for(i in 1:10000){
  # generate a candidate solution
  prop_x <- perturb_x(cur_x)
  # evaluate the candidate solution
  prop_fx <- evaluate_x(dists = y, tour = prop_x)
  # calculate the probability of accepting the candidate
  anneal_temp <- start_temp * temp_factor ^ i 
  accept_prob <- exp(-(prop_fx - cur_fx) / anneal_temp)
  # accept or reject the candidate
  if(prop_fx < cur_fx){ 
    cur_x <- prop_x 
    cur_fx <- prop_fx
  }
  else{ if(runif(1) < accept_prob){
    cur_x <- prop_x
    cur_fx <- prop_fx
  }
  }
  # store all results
  all_fx <- c(all_fx, cur_fx) 
  all_x <- rbind(all_x,cur_x)
  all_temp <- c(all_temp, anneal_temp)
}

# plot
pdf(file = "TSPCase12.pdf", width = 7, height = 7, compress = F)
plot(all_fx,type="l", xlab="Iterations", ylab="Distance", col = "darkred") 
dev.off()
# plot temp
plot(all_temp, type="l") # slower decrease 
# results 
best_solution <- which.min(all_fx) 
best_tour <- all_x[best_solution,] 
optimal_solution <- all_fx[best_solution] 
optimal_solution #12 842
best_tour
colnames(y)[as.numeric(best_tour)]

###############################
## Case 2.1
# use swap perturbation function
# swap operator
perturb_x2 <- function(tour){
  # select two cities at random
  ch <- sort(sample(2:(length(tour)-1),2,replace=FALSE)) 
  # swap them in the chain 
  p1 <- tour[ch[1]]
  tour[ch[1]] <- tour[ch[2]]
  tour[ch[2]] <- p1
  return(tour)
}
# set start temperature and geometric cooling factor 
set.seed(100) # for repeatability
start_temp <- 50000
temp_factor <- 0.999
# get an initial solution
cur_x <- get_initial_x(ncol(y))
# evaluate the solution
cur_fx <- evaluate_x(dists = y, tour = cur_x)
# initialize results data frames
all_fx <- c()
all_x <- data.frame()
all_temp <- c()

# for a fixed number of iterations
set.seed(2020)
for(i in 1:10000){
  # generate a candidate solution
  prop_x <- perturb_x2(cur_x)
  # evaluate the candidate solution
  prop_fx <- evaluate_x(dists = y, tour = prop_x)
  # calculate the probability of accepting the candidate
  anneal_temp <- start_temp * temp_factor ^ i 
  accept_prob <- exp(-(prop_fx - cur_fx) / anneal_temp)
  # accept or reject the candidate
  if(prop_fx < cur_fx){ 
    cur_x <- prop_x 
    cur_fx <- prop_fx
  }
  else{ if(runif(1) < accept_prob){
    cur_x <- prop_x
    cur_fx <- prop_fx
  }
  }
  # store all results
  all_fx <- c(all_fx, cur_fx) 
  all_x <- rbind(all_x,cur_x)
  all_temp <- c(all_temp, anneal_temp)
}

# plot
plot(all_fx,type="l") # convergence is quick 
plot(all_temp,type="l")
# results 
best_solution <- which.min(all_fx) 
best_tour <- all_x[best_solution,] 
optimal_solution <- all_fx[best_solution] 
optimal_solution # 13 363
best_tour
colnames(y)[as.numeric(best_tour)]
############################
## Case 3.1
# same settings as case 1.2 but log cooling schedule
set.seed(100) # for repeatability
start_temp <- 50000
temp_factor <- 0.999
# get an initial solution
cur_x <- get_initial_x(ncol(y))
# evaluate the solution
cur_fx <- evaluate_x(dists = y, tour = cur_x)
# initialize results data frames
all_fx <- c()
all_x <- data.frame()
all_temp <- c()

# for a fixed number of iterations
set.seed(2020)
for(i in 1:10000){
  # generate a candidate solution
  prop_x <- perturb_x(cur_x)
  # evaluate the candidate solution
  prop_fx <- evaluate_x(dists = y, tour = prop_x)
  # calculate the probability of accepting the candidate
  anneal_temp <- start_temp/(1+temp_factor*log(1+i))
  accept_prob <- exp(-(prop_fx - cur_fx) / anneal_temp)
  # accept or reject the candidate
  if(prop_fx < cur_fx){ 
    cur_x <- prop_x 
    cur_fx <- prop_fx
  }
  else{ if(runif(1) < accept_prob){
    cur_x <- prop_x
    cur_fx <- prop_fx
  }
  }
  # store all results
  all_fx <- c(all_fx, cur_fx) 
  all_x <- rbind(all_x,cur_x)
  all_temp <- c(all_temp, anneal_temp)
}

# plot
plot(all_fx,type="l") # convergence is quick 
plot(all_temp, type="l")
# results 
best_solution <- which.min(all_fx) 
best_tour <- all_x[best_solution,] 
optimal_solution <- all_fx[best_solution] 
optimal_solution # 20 955 - worse 
best_tour
colnames(y)[as.numeric(best_tour)]
############################
## Case 3.2
# same settings as case 1.2 but log cooling schedule
set.seed(100) # for repeatability
start_temp <- 50000
temp_factor <- 0.999
# get an initial solution
cur_x <- get_initial_x(ncol(y))
# evaluate the solution
cur_fx <- evaluate_x(dists = y, tour = cur_x)
# initialize results data frames
all_fx <- c()
all_x <- data.frame()
all_temp <- c()

# for a fixed number of iterations
set.seed(2020)
for(i in 1:10000){
  # generate a candidate solution
  prop_x <- perturb_x(cur_x)
  # evaluate the candidate solution
  prop_fx <- evaluate_x(dists = y, tour = prop_x)
  # calculate the probability of accepting the candidate
  anneal_temp <- start_temp/(1+temp_factor*(i^2))
  accept_prob <- exp(-(prop_fx - cur_fx) / anneal_temp)
  # accept or reject the candidate
  if(prop_fx < cur_fx){ 
    cur_x <- prop_x 
    cur_fx <- prop_fx
  }
  else{ if(runif(1) < accept_prob){
    cur_x <- prop_x
    cur_fx <- prop_fx
  }
  }
  # store all results
  all_fx <- c(all_fx, cur_fx) 
  all_x <- rbind(all_x,cur_x)
  all_temp <- c(all_temp, anneal_temp)
}

# plot
plot(all_fx,type="l") # convergence is quick 
plot(all_temp, type="l")
# results 
best_solution <- which.min(all_fx) 
best_tour <- all_x[best_solution,] 
optimal_solution <- all_fx[best_solution] 
optimal_solution # 12 984
best_tour
colnames(y)[as.numeric(best_tour)]
########## END 
