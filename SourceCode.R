# SIM AND OPT ASS1
rm(list = ls(all=TRUE))
setwd("~/Documents/Masters 2020:2021/Simulation and Optimisation/Optimisation/Ass1")
library(Rglpk)
library(xtable)
#######################
## Question 1
# demand
D <- c(1200,	1100,	1300,	1000)
# dvs
dv <- paste(c(rep("x", 8), rep("y", 12)),
            c(11,61,12,22,33,43,44,54,11,51,61,12,22,62,23,33,43,34,44,54),
            sep="")
# demand constr
demx <- as.matrix(simple_triplet_matrix(i = c(1,1,2,2,3,3,4,4), 
                                        j = c(1:8), 
                                        v = rep(1, 8), ncol = 20))
demy <- as.matrix(simple_triplet_matrix(i = c(1,1,1,2,2,2,3,3,3,4,4,4), 
                                        j = c(9:20), 
                                        v = rep(1, 12), ncol = 20))
dem <- rbind(demx, demy)
colnames(dem) <- dv

# dir
dem_sig <- rep(">=", nrow(dem))
# rhs
dem_rhs <- c(0.35*D, 0.65*D)

# manager's constr
# sumj(xij) - 3sumj(yij) <=0
manx <- as.matrix(simple_triplet_matrix(i = c(1,1,2,3,4,4,5,6), 
                                        j = c(1,3,4,5,6,7,8,2), 
                                        v = rep(1, 8), ncol = 8))
many <- as.matrix(simple_triplet_matrix(i = c(1,1,2,2,3,3,4,4,5,5,6,6), 
                                        j = c(1,4,5,7,8,10,9,11,2,12,3,6), 
                                        v = rep(-3, 12), ncol = 12))
# combine
man <-cbind(manx, many)
colnames(man) <- dv
# rhs
man_rhs <- rep(0, nrow(man))
# dir
man_sig <- rep("<=", nrow(man))

##########
# solve
coeff <- rbind(dem, man)
rhs <- c(dem_rhs, man_rhs)
dir <- c(dem_sig, man_sig)

# obj coeff
obj1 <- c(2,3,2,2.5,2,1,1,1.5,
          1.5,2,2.5,1.5,1,2.5,1,2,1,2,1,2)
names(obj1) <- dv

res1 <- Rglpk_solve_LP(obj = obj1, 
                       mat = coeff, 
                       rhs = rhs, 
                       dir = dir, 
                       max = FALSE)
# cost
res1$optimum
# sol
sol1 <- cbind(dv, res1$solution)[-which(res1$solution == 0),] 
colnames(sol1) <- c("dv", "value (tonnes)")
xtable(sol1)

###########################################
## Question 2
# new dv's: alphas, betas, d's
dv2 <- paste(c(rep("x", 8), rep("y", 12), rep("alpha", 6), rep("beta",6), rep("d",6)),
             c(11,61,12,22,33,43,44,54,11,51,61,12,22,62,23,33,43,34,44,54,rep(1:6,2), 1:6),
             sep="")
##### 
# demand constr
# fuel - x
demx <- as.matrix(simple_triplet_matrix(i = c(rep(1,2), rep(2,2), 
                                              rep(3,2), rep(4,2)), 
                                        j = c(1:8), 
                                        v = rep(1, 8), ncol = 38))
colnames(demx) <- dv2
# corn - y
demy <- as.matrix(simple_triplet_matrix(i = c(rep(1,3), rep(2,3), 
                                              rep(3,3), rep(4,3)), 
                                        j = c(9:20), 
                                        v = rep(1, 12), ncol = 38))
colnames(demy) <- dv2

demytot <- as.matrix(simple_triplet_matrix(i = c(rep(1,4), rep(2,4), rep(3,4), 
                                                 rep(4,4), rep(5,4), rep(6,4)), 
                                           j = c(9,12,21,27,
                                                 13,15,22,28,
                                                 16,18,23,29,
                                                 17,19,24,30,
                                                 10,20,25,31,
                                                 11,14,26,32), 
                                           v = rep(c(-1,-1,1,1),6), ncol = 38))
colnames(demytot) <- dv2
# combine all
dem <- rbind(demx, demy, demytot)
colnames(dem) <- dv2
# rhs
dem_rhs <- c(0.35*D, 0.65*D, rep(0, nrow(demytot)))
#signs
dem_sig <- c(rep(">=", 2*nrow(demx)), rep("==", nrow(demytot)))

# managers constr - amended
man2 <- cbind(man, matrix(0, ncol = 18, nrow=nrow(man)))
colnames(man2) <- dv2

# constr on alphas, betas (binary constrS)
alphs1 <- cbind(matrix(0,ncol=20, nrow = 6), diag(1, 6), 
                matrix(0,ncol=6, nrow = 6), diag(1000,6))
colnames(alphs1) <- dv2

bets1 <- cbind(matrix(0,ncol=26, nrow = 6), diag(1, 6), diag(-1000,6))
colnames(bets1) <- dv2

bets2 <- cbind(matrix(0, ncol=26, nrow=6), diag(1,6), diag(-10^6, 6))
colnames(bets2) <- dv2

# combine all binary constr
bin <- rbind(alphs1, bets1, bets2)
colnames(bin) <- dv2
# rhs
bin_rhs <- c(rep(1000, 6), rep(0,6), rep(0,6))
# dir
bin_sig <- c(rep("<=",6), rep(">=",6), rep("<=",6))
###########
# solve
coeff2 <- rbind(dem, man2, bin)
rhs2 <- c(dem_rhs, man_rhs, bin_rhs)
dir2 <- c(dem_sig, man_sig, bin_sig)

# obj coeff
ycost <- c(1.5,1,2,1,2,2.5)
obj2 <-c(2,3,2,2.5,2,1,1,1.5, # xij
         rep(0,12) , #yij
         ycost, #alphas
         0.75*ycost, #betas
         rep(0,6) #b's
)
names(obj2) <- dv2

res2 <- Rglpk_solve_LP(obj = obj2, 
                       mat = coeff2, 
                       rhs = rhs2, 
                       dir = dir2, 
                       max = FALSE, 
                       types = c(rep("C",32), rep("B",6)))
# cost
res2$optimum
sol2 <- cbind(dv2, res2$solution)[-which(res2$solution==0),]
colnames(sol2) <- c("dv", "value")
xtable(sol2)

#############################################
## Question 3
# dvs
dv3 <- paste(c(rep("x", 8), rep("y", 12), rep("g", 6), rep("x",6)),
             c(11,61,12,22,33,43,44,54,11,51,61,12,22,62,23,33,43,34,44,54,rep(1:6,2)),
             sep="")
# demand constr
demx <- as.matrix(simple_triplet_matrix(i = c(1,1,2,2,3,3,4,4), 
                                        j = c(1:8), 
                                        v = rep(1, 8), ncol = 32))
demy <- as.matrix(simple_triplet_matrix(i = c(1,1,1,2,2,2,3,3,3,4,4,4), 
                                        j = c(9:20), 
                                        v = rep(1, 12), ncol = 32))
dem <- rbind(demx, demy)
colnames(dem) <- dv3

# dir
dem_sig <- rep(">=", 8)
# rhs
dem_rhs <- c(0.35*D, 0.65*D)

# manager's constr
# sumj(xij) - 3sumj(yij) <=0
manx <- as.matrix(simple_triplet_matrix(i = c(1,1,2,2,3,3,4,4), 
                                        j = c(1:8), 
                                        v = rep(1,8), ncol = 8))
many <- as.matrix(simple_triplet_matrix(i = c(1,1,1,2,2,2,3,3,3,4,4,4), 
                                        j = c(1:12), 
                                        v = rep(-3,12), ncol = 12))
man <- cbind(manx, many, matrix(0, ncol=12, nrow = nrow(manx)))
colnames(man) <- dv3
# rhs
man_rhs <- rep(0, nrow(man))
# dir
man_sig <- rep("<=", nrow(man))

# total fuel constraint
# sum(xij) - xi = 0
totfuel <- as.matrix(simple_triplet_matrix(i = c(rep(1,3), rep(2,2), rep(3,2), 
                                                 rep(4,3), rep(5,2), rep(6,2)), 
                                           j = c(1,3,27,
                                                 4,28,
                                                 5,29,
                                                 6,7,30,
                                                 8,31,
                                                 2,32), 
                                           v = c(1,1,-1,
                                                 1,-1,
                                                 1,-1,
                                                 1,1,-1,
                                                 1,-1,
                                                 1,-1), ncol = 32))
colnames(totfuel) <- dv3
# rhs 
totfuel_rhs <- rep(0, nrow(totfuel))
# dir 
totfuel_sig <- rep("==", nrow(totfuel))

# binary constraint 
bin <- cbind(matrix(0, ncol = 20, nrow=6), 
             diag(-520, nrow = 6),
             diag(1, nrow = 6))
colnames(bin) <- dv3
# rhs
bin_rhs <- rep(0, nrow(bin))
# dir
bin_sig <- rep(">=", nrow(bin))

# previous month constraint
# xi - U*gammai-1 <= 0
U <- 10^7
prev <- cbind(matrix(0, ncol = 20, nrow=6), 
              cbind(rbind(rep(0,5),diag(U, 5)), c(U, rep(0,5))),
              diag(1, nrow = 6))
colnames(prev) <- dv3
# rhs
prev_rhs <- rep(U, nrow(prev))
# dir
prev_sig <- rep("<=", nrow(prev))
################
# solve
coeff <- rbind(dem, man, totfuel, bin)
rhs <- c(dem_rhs, man_rhs, totfuel_rhs, bin_rhs)
dir <- c(dem_sig, man_sig, totfuel_sig, bin_sig)

# obj coeff
obj3 <- c(2,3,2,2.5,2,1,1,1.5,
          1.5,2,2.5,1.5,1,2.5,1,2,1,2,1,2, 
          rep(0, 12))
names(obj3) <- dv3

res3 <- Rglpk_solve_LP(obj = obj3, 
                       mat = coeff, 
                       rhs = rhs, 
                       dir = dir, 
                       max = FALSE)
# cost
res3$optimum
sol3 <- cbind(dv3,res3$solution)[-which(res3$solution < 0.01),]
colnames(sol3) <- c("dv", "value")
xtable(sol3)
