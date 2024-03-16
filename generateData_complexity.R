#Modelle: 
#1 random intercept                       y = b0 + b1*obs + b2*cond + (1|subj) + epsilon
#1 random intercept + 1 random slope      y = b0 + b1*obs + b2*cond + b3 * group1 + (1 + group1|subj) + epsilon
#1 random intercepts + 2 random slopes    y = b0 + b1*obs + b2*cond + b3 * group1 + b4 * group2 + (1 + group1 + group2|subj) + epsilon
#nested (komplexer?)                      y = b0 + b1*obs + b2*cond + (1|subj) + (1|subj:group) + epsilon
#Interaktion?

library(MASS)
library(lme4)
library(faux)


sim_data_complexity <- function(n.obs = 1, n.pig = 24, n.test = 4, n.batch = 2, n.day = 4, n.pen = 4) {
  
  pig <- 1:n.pig
  test <- 1:n.test
  batch <- 1:n.batch
  day <- 1:n.day
  pen <- 1:n.pen
  sow <- 1:(n.pig/8)
  group <- 1:(n.pig/8)
  n.sow <- length(sow)
  n.group <- length(group)
  
  pig.d <- rep(pig, each = length(test) * n.obs)
  test.d <- rep(test, length(pig.d)/length(test))
  batch.d <- rep(batch, each = length(pig.d)/length(batch))
  day.d <- rep(rep(day, each = length(pig)*n.obs/length(day)), length(day))
  pen.d <- rep(rep(pen, each = length(test)*n.obs), length(pig.d) / length(test) / length(pen) / n.obs)
  sow.d <- rep(c(rep(sow, each = length(test) * n.obs ), rep(max(sow):1, each = length(test) * n.obs)), length(pig.d) / 2 / length(sow) / length(test) / n.obs)
  #group.d <- c(rep(rep(group, each = length(test) * n.obs), 2 * length(pig.d) / length(day) / length(group) / length(test) / n.obs), rep(rep(max(group):1, each = length(test) * n.obs), 2 * length(pig.d) / length(day) / length(group) / length(test) / n.obs))
  group.d <- rep(group, each = n.obs * n.test * n.pig / n.group)
  daytest.d <- as.factor(paste0(day.d, test.d))
  
  sd_pig <- 5
  sd_group <- 5
  sd_sow <- 5
  sd_pen <- 5
  
  int_pig <- rnorm(n.pig, 0, sd_pig)[pig.d]
  int_group <- rnorm(n.group, 0, sd_group)[group.d]
  int_sow <- rnorm(n.sow, 0, sd_sow)[sow.d]
  int_pen <- rnorm(n.pen, 0, sd_pen)[pen.d]
  
  slope_group.day_mat <- mvrnorm(n.group, mu = rep(0, n.day), Sigma = matrix(3, nrow = n.day, ncol = n.day) + diag(1, nrow = n.day))
  slope_group.day <- vector("numeric", length = length(pig.d))
  for(i in 1:length(pig.d)) {
    slope_group.day[i] <- slope_group.day_mat[group.d[i], day.d[i]]  
  }
  
  slope_group.test_mat <- mvrnorm(n.group, mu = rep(0, n.test), Sigma = matrix(3, nrow = n.test, ncol = n.test) + diag(1, nrow = n.test))
  slope_group.test <- vector("numeric", length = length(pig.d))
  for(i in 1:length(pig.d)) {
    slope_group.test[i] <- slope_group.test_mat[group.d[i], test.d[i]]  
  }
  
  slope_sow.test_mat <- mvrnorm(n.sow, mu = rep (0, n.test), Sigma = matrix(3, nrow = n.test, ncol = n.test) + diag(1, nrow = n.test))
  slope_sow.test <- vector("numeric", length = length(pig.d))
  for(i in 1:length(pig.d)) {
    slope_sow.test[i] <- slope_sow.test_mat[sow.d[i], test.d[i]]  
  }
  
  slope_sow.day_mat <- mvrnorm(n.sow, mu = rep (0, n.day), Sigma = matrix(3, nrow = n.day, ncol = n.day) + diag(1, nrow = n.day))
  slope_sow.day <- vector("numeric", length = length(pig.d))
  for(i in 1:length(pig.d)) {
    slope_sow.day[i] <- slope_sow.day_mat[sow.d[i], day.d[i]]  
  }
  
  slope_pen.daytest_mat <- mvrnorm(n.pen, mu = rep (0, n.day*n.test), Sigma = matrix(3, nrow = n.day*n.test, ncol = n.day*n.test) + diag(1, nrow = n.day*n.test))
  slope_pen.daytest <- vector("numeric", length = length(pig.d))
  for(i in 1:length(pig.d)) {
    slope_pen.daytest[i] <- slope_pen.daytest_mat[pen.d[i], daytest.d[i]]  
  }
  
  slope_pen.test_mat <- mvrnorm(n.pen, mu = rep (0, n.test), Sigma = matrix(3, nrow = n.test, ncol = n.test) + diag(1, nrow = n.test))
  slope_pen.test <- vector("numeric", length = length(pig.d))
  for(i in 1:length(pig.d)) {
    slope_pen.test[i] <- slope_pen.test_mat[pen.d[i], test.d[i]]  
  }
  
  slope_pen.day_mat <- mvrnorm(n.pen, mu = rep (0, n.day), Sigma = matrix(3, nrow = n.day, ncol = n.day) + diag(1, nrow = n.day))
  slope_pen.day <- vector("numeric")
  for(i in 1:length(pig.d)) {
    slope_pen.day[i] <- slope_pen.day_mat[pen.d[i], day.d[i]]  
  }
  
  slope_pen.batch_mat <- mvrnorm(n.pen, mu = rep (0, n.batch), Sigma = matrix(3, nrow = n.batch, ncol = n.batch) + diag(1, nrow = n.batch))
  slope_pen.batch <- vector("numeric", length = length(pig.d))
  for(i in 1:length(pig.d)) {
    slope_pen.batch[i] <- slope_pen.batch_mat[pen.d[i], batch.d[i]]  
  }
  
  return(data.frame(pig.d, test.d, batch.d, day.d, pen.d, sow.d, group.d, daytest.d,
                    int_pig, int_group, int_pen, int_sow,
                    slope_group.day, slope_group.test, slope_pen.batch, slope_pen.day, slope_pen.daytest, slope_pen.test,
                    slope_sow.day, slope_sow.test))
}


data <- read.csv("D:/Statistik BA/Assistenz Futschik/emmeans/Dataset_Exploring.csv")
table(data$pig)
data

pig <- 1:n.pig
test <- 1:n.test
batch <- 1:n.batch
day <- 1:n.day
pen <- 1:n.pen
sow <- 1:(n.pig/8)
group <- 1:(n.pig/8)
n.sow <- length(sow)
n.group <- length(group)

pig.d <- rep(pig, each = length(test) * n.obs)
test.d <- rep(test, length(pig.d)/length(test))
batch.d <- rep(batch, each = length(pig.d)/length(batch))
day.d <- rep(rep(day, each = length(pig)*n.obs/length(day)), length(day))
pen.d <- rep(rep(pen, each = length(test)*n.obs), length(pig.d) / length(test) / length(pen) / n.obs)
sow.d <- rep(c(rep(sow, each = length(test) * n.obs ), rep(max(sow):1, each = length(test) * n.obs)), length(pig.d) / 2 / length(sow) / length(test) / n.obs)
#group.d <- c(rep(rep(group, each = length(test) * n.obs), 2 * length(pig.d) / length(day) / length(group) / length(test) / n.obs), rep(rep(max(group):1, each = length(test) * n.obs), 2 * length(pig.d) / length(day) / length(group) / length(test) / n.obs))
group.d <- rep(group, each = n.obs * n.test * n.pig / n.group)
daytest.d <- as.factor(paste0(day.d, test.d))

sd_pig <- 5
sd_group <- 5
sd_sow <- 5
sd_pen <- 5

int_pig <- rnorm(n.pig, 0, sd_pig)[pig.d]
int_group <- rnorm(n.group, 0, sd_group)[group.d]
int_sow <- rnorm(n.sow, 0, sd_sow)[sow.d]
int_pen <- rnorm(n.pen, 0, sd_pen)[pen.d]

return(data.frame(pig.d, test.d, batch.d, day.d, pen.d, sow.d, group.d, daytest.d,
                  int_pig, int_group, int_pen, int_sow,
                  slope_group.day, slope_group.test, slope_pen.batch, slope_pen.day, slope_pen.daytest, slope_pen.test,
                  slope_sow.day, slope_sow.test))

sim_data_complexity <- function(n.obs = 5, n.ID = 60, n.x1 = 10, n.x2 = 6, n.x3 = 3) {
  ID <- 1:n.ID
  x1 <- 1:n.x1
  x2 <- 1:n.x2
  x3 <- 1:n.x3
  
  ID.d <- rep(ID, each = n.obs)
  x1.d <- rep(rep(x1, each = n.obs), length(ID.d) / n.obs / n.x1)
  x2.d <- rep(rep(x2, each = length(ID.d) / n.x2 / 4), 4)
  x3.d <- rep(x3, length(ID.d) / n.x3)
  data <- data.frame(ID = ID.d, x1 = x1.d, x2 = x2.d, x3 = x3.d)

  RE <- mvrnorm(n = n.ID, mu = rep(0, ncol(data)), Sigma = matrix(5, nrow = ncol(data), ncol = ncol(data)) + diag(10, nrow = ncol(data)))
  int_ID <- RE[ID.d, 1]
  slope_ID.x1 <- RE[ID.d, 2]
  slope_ID.x2 <- RE[ID.d, 3]
  slope_ID.x3 <- RE[ID.d, 4]
  
  data$y <- 10 + (2 + slope_ID.x1) * data$x1 + (5 + slope_ID.x3) * x3.d + int_ID + rnorm(length(data$ID), 0, 4)
  #return(data)
  return(cbind(data, int_ID, slope_ID.x1, slope_ID.x2, slope_ID.x3))
}

data <- sim_data_complexity(n.obs = 3, n.ID = 18, n.x1 = 3, n.x2 = 3, n.x3 = 3)

lmer(y ~ x1 + x3 + (1 + x1 + x3| ID), data = data)
sd(data$int_ID)
sd(data$slope_ID.x1)
max(cor(data[,1:4])[cor(data[,1:4]) != 1])


sim_data_complexity <- function(n.ID = 60, n.x1 = 10, n.x2 = 6, n.x3 = 5) {
  ID <- 1:n.ID
  x1 <- 1:n.x1
  x2 <- 1:n.x2
  x3 <- 1:n.x3
  
  ID.d <- rep(ID, each = n.x1 * n.x2 * n.x3)
  x1.d <- rep(rep(x1, each = n.x2 * n.x3), n.ID)
  x2.d <- rep(rep(x2, each = n.x1), n.x2 * n.x3)
  x3.d <- rep(rep(x3, n.x2), n.x1 * n.x2)
  data <- data.frame(ID = ID.d, x1 = x1.d, x2 = x2.d, x3 = x3.d)
  
  RE <- mvrnorm(n = n.ID, mu = rep(0, ncol(data)), Sigma = matrix(5, nrow = ncol(data), ncol = ncol(data)) + diag(30, nrow = ncol(data)))
  int_ID <- RE[ID.d, 1]
  slope_ID.x1 <- RE[ID.d, 2]
  slope_ID.x2 <- RE[ID.d, 3]
  slope_ID.x3 <- RE[ID.d, 4]
  
  data$y <- 10 + (2 + slope_ID.x1) * data$x1 + (4 + slope_ID.x2) * data$x2 + (5 + slope_ID.x3) * x3.d + int_ID + rnorm(length(data$ID), 0, 8)
  #return(data)
  return(cbind(data, int_ID, slope_ID.x1, slope_ID.x2, slope_ID.x3))
}

data <- sim_data_complexity(n.ID = 8, n.x1 = 8, n.x2 = 8, n.x3 = 8)
cor(data)
lmer(y ~ x1 + x2 + x3 + (1 + x1 + x2 + x3 | ID), data = data)
sd(data$slope_ID.x2)

