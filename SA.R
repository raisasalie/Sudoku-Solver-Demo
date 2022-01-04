
# SA 
rm(list = ls(all=T))
setwd("~/Documents/Masters 2020:2021/Simulation and Optimisation/Optimisation/Ass2")
# Solving Sudoku
# create a Sudoku puzzle
s <- matrix(0,ncol=9,nrow=9) 
s[1,c(6,8)] <- c(6,4) 
s[2,c(1:3,8)] <- c(2,7,9,5) 
s[3,c(2,4,9)] <- c(5,8,2) 
s[4,3:4] <- c(2,6) 
s[6,c(3,5,7:9)] <- c(1,9,6,7,3) 
s[7,c(1,3:4,7)] <- c(8,5,2,4) 
s[8,c(1,8:9)] <- c(3,8,5) 
s[9,c(1,7,9)] <- c(6,9,1)
s

# free spaces
free_spaces <- (s == 0) 
free_spaces

# get init sol
get_initial_x <- function(s){
  # identify free spaces
  free_spaces <- (s == 0)
  # fixed non-free spaces
  cur_x <- s
  # randomly choose a number between 1 and 9 for free spaces 
  cur_x[free_spaces] <- sample(1:9, sum(free_spaces), replace=T)
  return(cur_x) 
}

# init sol
init_x <- get_initial_x(s) 
init_x

# function to count duplicates
count_duplicates <- function(x) {
  n_dup <- length(x) - length(unique(x))
  return(n_dup) 
}

# eval solution
evaluate_x <- function(x){ # within-row duplications
  row_dups <- sum(apply(x,1,count_duplicates)) # within-col duplications
  col_dups <- sum(apply(x,2,count_duplicates))
  # within-block duplications
  block_dups <- 0 
  for (i in 1:3){
    for(j in 1:3){
      small_x <- x[(3*(i-1) + 1):(3*i), (3*(j-1) + 1):(3*j)] 
      thisblock_dups <- count_duplicates(as.vector(small_x)) 
      block_dups <- block_dups + thisblock_dups
    } 
  }
  total_dups <- row_dups + col_dups + block_dups
  return(total_dups)
}
# eval init sol
evaluate_x(init_x)

# perturb function
perturb_x <- function(x, free_spaces) {
  # select a free site
  site_to_change <- sample(1:81,1,prob=free_spaces) 
  # change that site at random
  x[site_to_change] <- sample(1:9,1,replace=T)
  return(x) 
}

# simulatiion
set.seed(100) # for repeatability 
start_temp <- 1e3
temp_factor <- 0.995

cur_x <- get_initial_x(s) 
cur_fx <- evaluate_x(cur_x)
cur_x
cur_fx

#############################
## Base Case 
# initial results data frames
all_fx <- c()
all_x <- data.frame()
# for a fixed number of iterations 
for(i in 1:3000){
# generate a candidate solution
prop_x <- perturb_x(cur_x, free_spaces = free_spaces) 
# evaluate the candidate solution
prop_fx <- evaluate_x(prop_x)
# calculate the probability of accepting the candidate
anneal_temp <- start_temp * temp_factor ^ i 
accept_prob <- exp(-(prop_fx - cur_fx) / anneal_temp)
# accept or reject the candidate
if(prop_fx < cur_fx){ 
    cur_x <- prop_x 
    cur_fx <- prop_fx
  }else{ 
    if(runif(1) < accept_prob){
    cur_x <- prop_x
    cur_fx <- prop_fx
    }}
# store all results
all_fx <- c(all_fx, cur_fx)
all_x <- rbind(all_x,as.vector(cur_x))
}

plot(all_fx,type="l")

plot(1:3000,start_temp * temp_factor ^ (1:3000))
tail(all_fx) # 23 errors
##########################################################
### CHANGES
#############################################
## Case 2.1
## logarithmic cooling sched
# initial results data frames
cur_x <- get_initial_x(s) 
all_fx <- c()
all_x <- data.frame()
# for a fixed number of iterations 
set.seed(100)
for(i in 1:3000){
  # generate a candidate solution
  prop_x <- perturb_x(cur_x, free_spaces = free_spaces) # evaluate the candidate solution
  prop_fx <- evaluate_x(prop_x)
  # calculate the probability of accepting the candidate
  anneal_temp <- start_temp/(1 +temp_factor*log(1+i))
  accept_prob <- exp(-(prop_fx - cur_fx) / anneal_temp)
  # accept or reject the candidate
  if(prop_fx < cur_fx){ 
    cur_x <- prop_x 
    cur_fx <- prop_fx
  }else{ 
    if(runif(1) < accept_prob){
      cur_x <- prop_x
      cur_fx <- prop_fx
    }}
  # store all results
  all_fx <- c(all_fx, cur_fx)
  all_x <- rbind(all_x,as.vector(cur_x))
}

plot(all_fx,type="l")

plot(1:3000,start_temp * temp_factor ^ (1:3000))
tail(all_fx) # 80 + errors

###########################################
## Case 2.2
## Quadratic cooling sched
# initial results data frames
cur_x <- get_initial_x(s) 
all_fx <- c()
all_x <- data.frame()
# for a fixed number of iterations 
set.seed(100)
for(i in 1:3000){
  # generate a candidate solution
  prop_x <- perturb_x(cur_x, free_spaces = free_spaces) # evaluate the candidate solution
  prop_fx <- evaluate_x(prop_x)
  # calculate the probability of accepting the candidate
  anneal_temp <- start_temp/(1+temp_factor*(i^2))
  accept_prob <- exp(-(prop_fx - cur_fx) / anneal_temp)
  # accept or reject the candidate
  if(prop_fx < cur_fx){ 
    cur_x <- prop_x 
    cur_fx <- prop_fx
  }else{ 
    if(runif(1) < accept_prob){
      cur_x <- prop_x
      cur_fx <- prop_fx
    }}
  # store all results
  all_fx <- c(all_fx, cur_fx)
  all_x <- rbind(all_x,as.vector(cur_x))
}

plot(all_fx,type="l")

plot(1:3000,start_temp * temp_factor ^ (1:3000))
tail(all_x)
cur_x
tail(all_fx) # 14 errors

########################################
## Case 3.1
# improve initial solution by row
get_initial_x_rows <- function(s){
  # identify free spaces
  free_spaces <- (s == 0)
  # initialize mnatrix to return
  cur_x <- matrix(NA, nrow = 9, ncol = 9)
  
  # for each row
  for (i in 1:nrow(s)){
    # the row
    row <- s[i,]
    # which spaces in row are free?
    free <- free_spaces[i,]
    # which numbers are not used in row?
    # if all zeroes , can use 1,..,9
    if(sum(free)==length(row)){not_taken <- c(1:9)}else{
      # if not all zeroes, can use those not in row
      not_taken <- c(1:9)[-row[!free]]
    }
    # replace open spaces with available numbers
    row[free] <- sample(x = not_taken, size = sum(free), replace = F)
    # replace in matrix to return 
    cur_x[i,] <- row
  }
  return(cur_x) 
}

# rerun with this adjustment
# simulatiion
set.seed(100) # for repeatability 
start_temp <- 1e3
temp_factor <- 0.995

cur_x <- get_initial_x_rows(s) 
cur_fx <- evaluate_x(cur_x)

# initial results data frames
all_fx <- c()
all_x <- data.frame()
# for a fixed number of iterations 
for(i in 1:6000){
  # generate a candidate solution
  prop_x <- perturb_x(cur_x, free_spaces = free_spaces) 
  # evaluate the candidate solution
  prop_fx <- evaluate_x(prop_x)
  # calculate the probability of accepting the candidate
  anneal_temp <- start_temp * temp_factor ^ i 
  accept_prob <- exp(-(prop_fx - cur_fx) / anneal_temp)
  # accept or reject the candidate
  if(prop_fx < cur_fx){ 
    cur_x <- prop_x 
    cur_fx <- prop_fx
  }else{ 
    if(runif(1) < accept_prob){
      cur_x <- prop_x
      cur_fx <- prop_fx
    }}
  # store all results
  all_fx <- c(all_fx, cur_fx)
  all_x <- rbind(all_x,as.vector(cur_x))
}

plot(all_fx,type="l")

plot(1:6000,start_temp * temp_factor ^ (1:6000))
tail(all_fx) # 12 errors - improvement

#######################################
## Case 3.2
# improve initial solution by block
get_initial_x_blocks <- function(s){
  # initialize mnatrix to return
  cur_x <- matrix(NA, nrow = 9, ncol = 9)
  
  for (b in 1:3){
    for (a in 1:3){
      # for each block
      block <- s[(3*b-2):(3*b),(3*a-2):(3*a)]
      # identify free spaces
      free <- (block == 0)
      # which numbers are not taken 
      # if all zeroes , can use 1,..,9
      if(sum(free)==9){not_taken <- c(1:9)}else{
        # if not all zeroes, can use those not in row
        not_taken <- c(1:9)[-block[!free]]
      }
      # replace open spaces with available numbers
      block[free] <- sample(x = not_taken, size = sum(free), replace = F)
      # replace in matrix to return 
      cur_x[(3*b-2):(3*b),(3*a-2):(3*a)] <- block
    }
  }
  return(cur_x) 
}

# rerun with this adjustment
# simulatiion
set.seed(100) # for repeatability 
start_temp <- 1e3
temp_factor <- 0.995

cur_x <- get_initial_x_blocks(s) 
cur_fx <- evaluate_x(cur_x)

# initial results data frames
all_fx <- c()
all_x <- data.frame()
# for a fixed number of iterations 
for(i in 1:6000){
  # generate a candidate solution
  prop_x <- perturb_x(cur_x, free_spaces = free_spaces) 
  # evaluate the candidate solution
  prop_fx <- evaluate_x(prop_x)
  # calculate the probability of accepting the candidate
  anneal_temp <- start_temp * temp_factor ^ i 
  accept_prob <- exp(-(prop_fx - cur_fx) / anneal_temp)
  # accept or reject the candidate
  if(prop_fx < cur_fx){ 
    cur_x <- prop_x 
    cur_fx <- prop_fx
  }else{ 
    if(runif(1) < accept_prob){
      cur_x <- prop_x
      cur_fx <- prop_fx
    }}
  # store all results
  all_fx <- c(all_fx, cur_fx)
  all_x <- rbind(all_x,as.vector(cur_x))
}

plot(all_fx,type="l")

plot(1:6000,start_temp * temp_factor ^ (1:6000))
tail(all_fx) # 10 errors

#######################################
## Case 4.1
# perturb function - second version
perturb_x2 <- function(x, free_spaces) {
  # select a free site
  # choose a row
  i <- sample(1:9,1)
  # choose a col (where there is a free site ij)
  j <- sample(1:9,1,prob=free_spaces[i,])
  # row of site
  row <- s[i,]
  # column of site
  col <- s[,j]
  # which numbers can we choose from?
  taken <- c(row, col)[-which(c(row,col)==0)]
  
  # change that site to a legal value given fixed sites
  x[i,j] <- sample(c(1:9)[-taken],1)
  
  return(x) 
}

# simulatiion
set.seed(100) # for repeatability 
start_temp <- 1e3
temp_factor <- 0.995

cur_x <- get_initial_x_blocks(s) 
cur_fx <- evaluate_x(cur_x)

# initial results data frames
all_fx <- c()
all_x <- data.frame()
# for a fixed number of iterations 
for(i in 1:3000){
  # generate a candidate solution
  prop_x <- perturb_x2(cur_x, free_spaces = free_spaces) 
  # evaluate the candidate solution
  prop_fx <- evaluate_x(prop_x)
  # calculate the probability of accepting the candidate
  anneal_temp <- start_temp * temp_factor ^ i 
  accept_prob <- exp(-(prop_fx - cur_fx) / anneal_temp)
  # accept or reject the candidate
  if(prop_fx < cur_fx){ 
    cur_x <- prop_x 
    cur_fx <- prop_fx
  }else{ 
    if(runif(1) < accept_prob){
      cur_x <- prop_x
      cur_fx <- prop_fx
    }}
  # store all results
  all_fx <- c(all_fx, cur_fx)
  all_x <- rbind(all_x,as.vector(cur_x))
}

plot(all_fx,type="l")

plot(1:3000,start_temp * temp_factor ^ (1:3000))
tail(all_fx) # 10 errors    
###############################
## Case 4.2
# as 4.1 with quadratic cooling sched 
# simulatiion
set.seed(100) # for repeatability 
start_temp <- 1e3
temp_factor <- 0.995

cur_x <- get_initial_x_blocks(s) 
cur_fx <- evaluate_x(cur_x)

# initial results data frames
all_fx <- c()
all_x <- data.frame()
# for a fixed number of iterations 
for(i in 1:6000){
  # generate a candidate solution
  prop_x <- perturb_x2(cur_x, free_spaces = free_spaces) 
  # evaluate the candidate solution
  prop_fx <- evaluate_x(prop_x)
  # calculate the probability of accepting the candidate
  anneal_temp <- start_temp/(1+temp_factor*(i^2))
  accept_prob <- exp(-(prop_fx - cur_fx) / anneal_temp)
  # accept or reject the candidate
  if(prop_fx < cur_fx){ 
    cur_x <- prop_x 
    cur_fx <- prop_fx
  }else{ 
    if(runif(1) < accept_prob){
      cur_x <- prop_x
      cur_fx <- prop_fx
    }}
  # store all results
  all_fx <- c(all_fx, cur_fx)
  all_x <- rbind(all_x,as.vector(cur_x))
}

plot(all_fx,type="l")

plot(1:6000,start_temp * temp_factor ^ (1:6000))
tail(all_fx) # 8 errors    
#########################################
## Case 4.3
# as 4.1 with quadratic cooling sched 
# try varying alphas
alphas <- seq(0.8, 0.999, length.out = 5)
# storage of results
res <- list()
for (k in 1:length(alphas)){
  # simulatiion
  set.seed(100) # for repeatability 
  start_temp <- 5e3
  temp_factor <- alphas[k]
  
  cur_x <- get_initial_x_blocks(s) 
  cur_fx <- evaluate_x(cur_x)
  
  # initial results data frames
  all_fx <- c()
  all_x <- data.frame()
  # for a fixed number of iterations 
  for(i in 1:8000){
    # generate a candidate solution
    prop_x <- perturb_x2(cur_x, free_spaces = free_spaces) 
    # evaluate the candidate solution
    prop_fx <- evaluate_x(prop_x)
    # calculate the probability of accepting the candidate
    anneal_temp <- start_temp/(1+temp_factor*(i^2))
    accept_prob <- exp(-(prop_fx - cur_fx) / anneal_temp)
    # accept or reject the candidate
    if(prop_fx < cur_fx){ 
      cur_x <- prop_x 
      cur_fx <- prop_fx
    }else{ 
      if(runif(1) < accept_prob){
        cur_x <- prop_x
        cur_fx <- prop_fx
      }}
    # store all results
    all_fx <- c(all_fx, cur_fx)
    all_x <- rbind(all_x,as.vector(cur_x))
  }
  # add to results
  res[[k]] <- all_fx
}

cols <- c("darkred", "red", "orange", "yellow", "green ",
          "blue", "darkblue", "purple", "violet", "pink")
pdf(file = "Case43.pdf", compress = F, width = 7, height=7)
plot(x = c(1:8000), y = res[[1]],type="l", col=cols[1], ylim = c(0,60), lwd=2, 
     xlab="Iteration", ylab="Errors")
for (i in 2:length(res)){
  lines(res[[i]], col=cols[i], lwd=2)
}
legend("topright",legend = paste(alphas), col = cols, lty = 1, lwd = 2)
dev.off()

# plot temps 
X <- c(1:5000)
plot(x = X, y = start_temp/(1+alphas[1]*(X^2)), type="l", col=cols[1], lwd=2,
     xlab="Run", ylab="Temperature", ylim=c(0,500))
for (i in 2:length(alphas)){
  lines(x=X, y = start_temp/(1+alphas[i]*(X^2)), lwd=2, col=cols[i])
}
#########################################
## Case 4.4
# choose best alpha = 0.94925

  # simulatiion
  set.seed(100) # for repeatability 
  start_temp <- 5e3
  temp_factor <- alphas[4]
  
  cur_x <- get_initial_x_blocks(s) 
  cur_fx <- evaluate_x(cur_x)
  
  # initial results data frames
  all_fx <- c()
  all_x <- data.frame()
  # for a fixed number of iterations 
  for(i in 1:15000){
    # generate a candidate solution
    prop_x <- perturb_x2(cur_x, free_spaces = free_spaces) 
    # evaluate the candidate solution
    prop_fx <- evaluate_x(prop_x)
    # calculate the probability of accepting the candidate
    anneal_temp <- start_temp/(1+temp_factor*(i^2))
    accept_prob <- exp(-(prop_fx - cur_fx) / anneal_temp)
    # accept or reject the candidate
    if(prop_fx < cur_fx){ 
      cur_x <- prop_x 
      cur_fx <- prop_fx
    }else{ 
      if(runif(1) < accept_prob){
        cur_x <- prop_x
        cur_fx <- prop_fx
      }}
    # store all results
    all_fx <- c(all_fx, cur_fx)
    all_x <- rbind(all_x,as.vector(cur_x))
  }

pdf(file = "Case44.pdf", compress = F, width = 7, height=7)
plot(x = c(1:15000), y = all_fx,type="l", col=cols[1], ylim = c(0,60), lwd=2, 
     xlab="Iteration", ylab="Errors")
dev.off()
all_fx[which.min(all_fx)]

###################################
## Case 5
# brute force
# simulatiion
set.seed(100) # for repeatability 
start_temp <- 1e3
temp_factor <- 0.995

cur_x <- get_initial_x(s) 
cur_fx <- evaluate_x(cur_x)

# initial results data frames
all_fx <- c()
all_x <- data.frame()
# for a fixed number of iterations 
for(i in 1:20000){
  # generate a candidate solution
  prop_x <- perturb_x2(cur_x, free_spaces = free_spaces) 
  # evaluate the candidate solution
  prop_fx <- evaluate_x(prop_x)
  # calculate the probability of accepting the candidate
  anneal_temp <- start_temp * temp_factor ^ i 
  accept_prob <- exp(-(prop_fx - cur_fx) / anneal_temp)
  # accept or reject the candidate
  if(prop_fx < cur_fx){ 
    cur_x <- prop_x 
    cur_fx <- prop_fx
  }else{ 
    if(runif(1) < accept_prob){
      cur_x <- prop_x
      cur_fx <- prop_fx
    }}
  # store all results
  all_fx <- c(all_fx, cur_fx)
  all_x <- rbind(all_x,as.vector(cur_x))
}

pdf(file = "Case5.pdf", width = 7, height = 7, compress = F)
plot(all_fx,type="l", xlab="Iteration", ylab = "Errors", lwd=2, col = cols[1])
dev.off()

plot(1:20000,start_temp * temp_factor ^ (1:20000))

all_fx[which.min(all_fx)]
# 4 errors with 20000 runs 
#save.image(file = "sudoku_pert2.RData")

#########################
## Case 6 - not included
## brute force with adjustments
# simulatiion
set.seed(100) # for repeatability 
start_temp <- 5e3
temp_factor <- alphas[4]

cur_x <- get_initial_x_blocks(s) 
cur_fx <- evaluate_x(cur_x)

# initial results data frames
all_fx2 <- c()
all_x2 <- data.frame()
# for a fixed number of iterations 
for(i in 1:20000){
  # generate a candidate solution
  prop_x <- perturb_x2(cur_x, free_spaces = free_spaces) 
  # evaluate the candidate solution
  prop_fx <- evaluate_x(prop_x)
  # calculate the probability of accepting the candidate
  anneal_temp <- start_temp/(1+temp_factor*(i^2))
  accept_prob <- exp(-(prop_fx - cur_fx) / anneal_temp)
  # accept or reject the candidate
  if(prop_fx < cur_fx){ 
    cur_x <- prop_x 
    cur_fx <- prop_fx
  }else{ 
    if(runif(1) < accept_prob){
      cur_x <- prop_x
      cur_fx <- prop_fx
    }}
  # store all results
  all_fx2 <- c(all_fx2, cur_fx)
  all_x2 <- rbind(all_x2, as.vector(cur_x))
}

pdf(file = "Case6.pdf", compress = F, width = 7, height=7)
plot(x = c(1:20000), y = all_fx2,type="l", col=cols[1], ylim = c(0,60), lwd=2, 
     xlab="Iteration", ylab="Errors")
dev.off()
all_fx2[which.min(all_fx2)]





