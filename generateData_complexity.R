#Modelle: 
#1 random intercept                       y = b0 + b1*obs + b2*cond + (1|subj) + epsilon
#1 random intercept + 1 random slope      y = b0 + b1*obs + b2*cond + b3 * group1 + (1 + group1|subj) + epsilon
#1 random intercepts + 2 random slopes    y = b0 + b1*obs + b2*cond + b3 * group1 + b4 * group2 + (1 + group1 + group2|subj) + epsilon
#nested (komplexer?)                      y = b0 + b1*obs + b2*cond + (1|subj) + (1|subj:group) + epsilon

#generatign data for models of different complexity
#n.schools ... number of school
#n.classes ... number of classes per school
#n.students ... number of students per class
# each student is either in condition 0 or 1
# all amounts must be even
#model:
##1 ... y = b0 + b_age * age + (1|obs)
##2 ... y = b0 + b_age * age + b_school * school + (1 + school|obs)
sim_data <- function(model, n.schools = 4, n.classes = 4, n.students = 20, n.obs = 4, b0 = 10, beta_age = 0, sd.int_obs = 5, sd_eps = 3) {
  student <- rep(1:(n.students*n.classes*n.schools), each = n.obs)
  observation <- as.factor(rep(1:n.obs, length(student)/n.obs))
  age <- sample(6:10, length(student), replace = TRUE)

  if(model == 1) {
    int_obs <- rep(rnorm(n.obs, 0, sd.int_obs), length(student)/n.obs)
    y <- b0 + beta_age * age + int_obs + rnorm(length(student), 0, sd_eps)
    return(data.frame(student, observation, age, y))
  }
  else if(model == 2) {
    school <- as.factor(rep(1:n.schools, each = length(student)/n.schools))
    class <- rep(rep(1:n.classes, each = n.obs), length(student)/(n.obs*n.classes))
    condition <- student %% 2
  }
}

data <- sim_data(model = 1, n.schools = 4, n.classes = 4, n.students = 10, n.obs = 6)
summary(lmer(y ~ age + (1|observation), data = data))
