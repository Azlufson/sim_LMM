#Modelle: 
#1 random intercept                       y = b0 + b1*obs + b2*cond + (1|subj) + epsilon
#1 random intercept + 1 random slope      y = b0 + b1*obs + b2*cond + b3 * group1 + (1 + group1|subj) + epsilon
#1 random intercepts + 2 random slopes    y = b0 + b1*obs + b2*cond + b3 * group1 + b4 * group2 + (1 + group1 + group2|subj) + epsilon
#nested (komplexer?)                      y = b0 + b1*obs + b2*cond + (1|subj) + (1|subj:group) + epsilon
#Interaktion?

library(MASS)
library(lme4)

#generatign data for models of different complexity
#n.schools ... number of school
#n.classes ... number of classes per school
#n.students ... number of students per class
# each student is either in condition 0 or 1
# all amounts must be even
#model:
##1 ... y = b0 + b_age * age + (1|obs)
##2 ... y = b0 + b_age * age + b_school * school + (1 + school|obs)
##3 ... y = b0 + b_age * age + b_school * school + b_class * class + (1 + school + class|obs)
sim_data <- function(model, n.schools = 4, n.classes = 4, n.students = 20, n.obs = 4, 
                     b0 = 10, beta_age = 0, sd.int_obs = 2, sd_eps = 1,
                     beta_school = 20, sd.slope_school = 10, corr_school.obs = .3,
                     corr_class.obs = .5, sd.slope_class = 15, beta_class = 5, corr_class.school = .8) {
  student <- rep(1:(n.students*n.classes*n.schools), each = n.obs)
  observation <- rep(1:n.obs, length(student)/n.obs)
  age <- sample(6:10, length(student), replace = TRUE)
  #condition <- student %% 2

  if(model == 1) {
    int_obs <- rep(rnorm(n.obs, 0, sd.int_obs), length(student)/n.obs)
    y <- b0 + beta_age * age + int_obs + rnorm(length(student), 0, sd_eps)
    return(data.frame(student, observation, age, y, int_obs))
  }
  else if(model == 2) {
    school <- rep(1:n.schools, each = length(student)/n.schools)

    cov_school.obs <- corr_school.obs * sd.int_obs * sd.slope_school #aus Korrelation und SDs Kovarianz der random effects berechnen
    sigma <- matrix(c(sd.int_obs^2, cov_school.obs,
                      cov_school.obs, sd.slope_school^2), nrow = 2, byrow = TRUE) #Kovarianzmatrix der random effects
    
    re_school.obs <- mvrnorm(n.obs, c(0,0), sigma)
    int_obs <- rep(re_school.obs[,1], n.students)
    slope_school <- rep(re_school.obs[,2], n.students)
    
    y <- b0 + beta_age * age + (beta_school + slope_school) * school + int_obs + rnorm(length(student), 0, sd_eps) #daten simulieren
    return(data.frame(student, observation, school, age, y, int_obs, slope_school))
  }
  
  else if (model == 3) {
    class <- rep(rep(1:n.classes, each = n.obs), length(student)/(n.obs*n.classes))
    school <- rep(1:n.schools, each = length(student)/n.schools)
    
    cov_class.obs <- corr_class.obs * sd.int_obs * sd.slope_class
    cov_school.obs <- corr_school.obs * sd.int_obs * sd.slope_school #aus Korrelation und SDs Kovarianz der random effects berechnen
    cov_class.school <- corr_class.school * sd.slope_school * sd.slope_class
    sigma <- matrix(c(sd.int_obs^2, cov_school.obs, cov_class.obs,
                      cov_school.obs, sd.slope_school^2, cov_class.school,
                      cov_class.obs, cov_class.school, sd.slope_class^2), nrow = 3, byrow = TRUE) #Kovarianzmatrix der random effects
    
    re_school.obs <- mvrnorm(n.obs, c(0,0, 0), sigma)
    int_obs <- rep(re_school.obs[,1], n.students)
    slope_school <- rep(re_school.obs[,2], n.students)
    slope_class <- rep(re_school.obs[,3], n.students)
    y <- b0 + beta_age * age + (beta_school + slope_school) * school + (beta_class + slope_class) * class + int_obs + rnorm(length(student), 0, sd_eps) #daten simulieren
    return(data.frame(student, observation, school, class, age, y, int_obs, slope_school))
    
  }
}


data <- read.csv("D:/Statistik BA/Assistenz Futschik/emmeans/Dataset_Exploring.csv")
table(data$pig)
data

lmm_full <- lmer(exploring ~ test*day + batch + contactbyhuman + (1|pig) + (1+day+test+contactbyhuman||sow) + (1+test+day+contactbyhuman||group) + (1+batch+day+test+test*day+contactbyhuman||pen), data = data)
lmm_full

pig <- 1:32
test <- 1:2
batch <- 1:2
day <- 1:2
pen <- 1:2
sow <- 1:(length(pig)/2)
group <- 1:(length(pig)/4)

pig.d <- rep(pig, each = length(test))
test.d <- rep(test, length(pig.d)/length(test))
batch.d <- rep(batch, each = length(pig.d)/length(batch))
day.d <- rep(rep(day, each = length(pig)/length(day)), length(day))
pen.d <- rep(rep(pen, each = length(test)), length(pig.d) / length(test) / length(pen))
sow.d <- rep(c(rep(sow, each = length(test)), rep(max(sow):1, each = length(test))), length(pig.d) / 2 / length(sow) / length(test))
group.d <- c(rep(rep(group, each = length(test)), length(pig.d) / length(day) / length(group) / length(test)), rep(rep(max(group):1, each = length(test)), length(pig.d) / length(day) / length(group) / length(test)))
daytest.d <- as.factor(paste0(day.d, test.d))

data_new <- data.frame(pig.d, test.d, batch.d, day.d, pen.d, sow.d, group.d, daytest.d)
data_new
table(data_new$sow.d)

sd_pig <- 1
sd_group <- 2
sd_sow <- 3
sd_pen <- 4
sd_testA_group <- 5
sd_testB_group <- 6
corr_testA.testB_group <- .1
cov_testA.testB_group <- corr_testA.testB_group * sd_testA_group * sd_testB_group
sd_day1_group <- 7
sd_day2_group <- 8
corr_day1.day2_group <- .2
cov_day1.day2_group <- corr_day1.day2_group * sd_day1_group * sd_day2_group
sd_testA_sow <- 9
sd_testB_sow <- 10
corr_testA.testB_sow <- .3
cov_testA.testB_sow <- corr_testA.testB_sow * sd_testA_sow * sd_testB_sow
sd_day1_sow <- 11
sd_day2_sow <- 12
cor_day1.day2_sow <- .4
cov_day1.day2_sow <- cor_day1.day2_sow * sd_day1_sow * sd_day2_sow
sd_day1testA_pen <- 13
sd_day2testA_pen <- 14
sd_day1testB_pen <- 15
sd_day2testB_pen <- 16
cor_day1testA.day2testA_pen <- .5
cor_day1testA.day1testB_pen <- .6
cor_day1testA.day2testB_pen <- .7
cor_day2testA.day1testB_pen <- .8
cor_day2testA.day2testB_pen <- .9
cor_day1testB.day2testB_pen <- -.1
cov_day1testA.day2testA_pen <- cor_day1testA.day2testA_pen * sd_day1testA_pen * sd_day1testB_pen
cov_day1testA.day1testB_pen <- cor_day1testA.day1testB_pen * sd_day1testA_pen * sd_day1testB_pen
cov_day1testA.day2testB_pen <- cor_day1testA.day2testB_pen * sd_day1testA_pen * sd_day2testB_pen
cov_day2testA.day1testB_pen <- cor_day2testA.day1testB_pen * sd_day2testA_pen * sd_day1testB_pen
cov_day2testA.day2testB_pen <- cor_day2testA.day2testB_pen * sd_day2testA_pen * sd_day2testB_pen
cov_day1testB.day2testB_pen <- cor_day1testB.day2testB_pen * sd_day1testB_pen * sd_day2testB_pen
sigma_daytest_pen <- matrix(c(sd_day1testA_pen^2, cor_day1testA.day2testA_pen, cor_day1testA.day1testB_pen, cor_day1testA.day2testB_pen, 
                              cor_day1testA.day2testA_pen, sd_day2testA_pen^2, cor_day2testA.day1testB_pen, cor_day2testA.day2testB_pen,
                              cor_day1testA.day1testB_pen, cor_day2testA.day1testB_pen, sd_day1testB_pen^2, cor_day1testB.day2testB_pen,
                              cor_day1testA.day2testB_pen, cor_day2testA.day2testB_pen, cor_day1testB.day2testB_pen, sd_day2testB_pen^2),
                            byrow = TRUE, nrow = 4)
sd_testA_pen <- 17
sd_testB_pen <- 18
cor_testA.testB_pen <- -.2
cov_testA.testB_pen <- cor_testA.testB_pen * sd_testA_pen * sd_testB_pen
sd_day1_pen <- 19
sd_day2_pen <- 20
cor_day1.day2_pen <- -.3
cov_day1.day2_pen <- cor_day1.day2_pen * sd_day1_pen * sd_day2_pen
sd_batch1_pen <- 21
sd_batch2_pen <- 22
cor_batch1.batch2_pen <- -.4
cov_batch1.batch2_pen <- cor_batch1.batch2_pen * sd_batch1_pen * sd_batch2_pen

int_pig <- rnorm(length(pig), 0, sd_pig)[pig.d]
int_group <- rnorm(length(group), 0, sd_group)[group.d]
int_sow <- rnorm(length(sow), 0, sd_sow)[sow.d]
int_pen <- rnorm(length(pen), 0, sd_pen)[pen.d]

slope_group.day_mat <- mvrnorm(length(group.d), mu = rep(0, length(day)), Sigma = matrix(c(sd_day1_group^2, cov_day1.day2_group,
                                                                                           cov_day1.day2_group, sd_day2_group^2), byrow = TRUE, nrow = 2))
slope_group.day <- vector("numeric")
for(i in 1:length(day.d)) {
  slope_group.day[i] <- slope_group.day_mat[group.d[i], day.d[i]]  
}

slope_group.test_mat <- mvrnorm(length(group), mu = rep(0, length(test)), Sigma = matrix(c(sd_testA_group^2, cov_testA.testB_group,
                                                                                   cov_testA.testB_group, sd_testB_group^2), byrow = TRUE, nrow = 2))
slope_group.test <- vector("numeric")
for(i in 1:length(day.d)) {
  slope_group.test[i] <- slope_group.test_mat[group.d[i], test.d[i]]  
}

slope_sow.test_mat <- mvrnorm(length(sow), mu = rep (0, length(test)), Sigma = matrix(c(sd_testA_sow^2, cov_testA.testB_sow,
                                                                                   cov_testA.testB_sow, sd_testB_sow^2), byrow = TRUE, nrow = 2))
slope_sow.test <- vector("numeric")
for(i in 1:length(day.d)) {
  slope_sow.test[i] <- slope_sow.test_mat[sow.d[i], test.d[i]]  
}

slope_sow.day_mat <- mvrnorm(length(sow), mu = rep (0, length(day)), Sigma = matrix(c(sd_day1_sow^2, cov_day1.day2_sow,
                                                                                        cov_day1.day2_sow, sd_day2_sow^2), byrow = TRUE, nrow = 2))
slope_sow.day <- vector("numeric")
for(i in 1:length(day.d)) {
  slope_sow.day[i] <- slope_sow.day_mat[sow.d[i], day.d[i]]  
}

slope_pen.daytest_mat <- mvrnorm(length(pen), mu = rep (0, length(day)*length(test)), Sigma = sigma_daytest_pen)
slope_pen.daytest <- vector("numeric")
for(i in 1:length(day.d)) {
  slope_pen.daytest[i] <- slope_pen.daytest_mat[pen.d[i], daytest[i]]  
}

slope_pen.test_mat <- mvrnorm(length(pen), mu = rep (0, length(test)), Sigma = matrix(c(sd_testA_pen^2, cov_testA.testB_pen,
                                                                                        cov_testA.testB_pen, sd_testB_pen^2), byrow = TRUE, nrow = 2))
slope_pen.test <- vector("numeric")
for(i in 1:length(day.d)) {
  slope_pen.test[i] <- slope_pen.test_mat[pen.d[i], test.d[i]]  
}

slope_pen.day_mat <- mvrnorm(length(pen), mu = rep (0, length(day)), Sigma = matrix(c(sd_day1_pen^2, cov_day1.day2_pen,
                                                                                      cov_day1.day2_pen, sd_day2_pen^2), byrow = TRUE, nrow = 2))
slope_pen.day <- vector("numeric")
for(i in 1:length(day.d)) {
  slope_pen.day[i] <- slope_pen.day_mat[pen.d[i], day.d[i]]  
}
cbind(pen.d, day.d, slope_pen.day)

slope_pen.batch_mat <- mvrnorm(length(pen), mu = rep (0, length(batch)), Sigma = matrix(c(sd_batch1_pen^2, cov_batch1.batch2_pen,
                                                                                      cov_batch1.batch2_pen, sd_batch2_pen^2), byrow = TRUE, nrow = 2))
slope_pen.batch <- vector("numeric")
for(i in 1:length(day.d)) {
  slope_pen.batch[i] <- slope_pen.batch_mat[pen.d[i], batch.d[i]]  
}


#exploring ~ test*day + batch + contactbyhuman + (1|pig) + (1+day+test+contactbyhuman||sow) + (1+test+day+contactbyhuman||group) + (1+batch+day+test+test*day+contactbyhuman||pen)
exploring <- 10 + 2 * test.d + 3 * day.d + 4 * as.numeric(daytest) + 5 * batch.d + int_pig + int_sow + (slope_sow.day + slope_sow.test) * sow.d + int_group + (slope_group.test + slope_group.day) * group + int_pen + (slope_pen.batch + slope_pen.day + slope_pen.test + slope_pen.daytest)



