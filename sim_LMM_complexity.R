#Alpha-Fehler, Güte:
#Modell erstellen
#testen: fixed effects, random effects
#Parameter: Stichprobengröße, Stärke des Effekts, Komplexität des Modells, missing values, balanciertes Design
#Schätzmethoden:  ML, REML
#                 Sattherwaite, Kenward-Rogers
#                 MCMC (Baayen et al., 2008)
#                 t as z

#mcmc nicht mehr unterstützt von lme4
#library(languageR)
#lmm <- lmer(model, sim_data_int())
#pvals.fnc(lmm)

####Modellkomplexität

##TODO: KR für REML?, PB für REML?
##      funktionen aus anderem skript importieren?
##      models für ungerade zahlen
##      future_sapply vs future_replicate testen
##      nicht nur factorial design?

#Modelle: 
#1 random intercept                       y = b0 + b1*obs + b2*cond + (1|subj) + epsilon
#1 random intercept + 1 random slope      y = b0 + b1*obs + b2*cond + b3 * group1 + (1 + group1|subj) + epsilon
#1 random intercepts + 2 random slopes    y = b0 + b1*obs + b2*cond + b3 * group1 + b4 * group2 + (1 + group1 + group2|subj) + epsilon
#nested (komplexer?)                      y = b0 + b1*obs + b2*cond + (1|subj) + (1|subj:group) + epsilon


library(future.apply)
library(parallel)
library(tidyverse)
library(lme4)
library(lmerTest)
library(pbkrtest)
library(afex)
library(MASS)

##Datengeneration
#einfaches Modell nur mit random intercept
# y = b0 + b1*obs + b2*cond + (1|subj) + epsilon
#n.subj und n.obs müssen gerade sein
# sim_data_1RI <- function(n.subj = 10, n.obs = 6, b0 = 10, beta_age = 0, beta_cond = 5, sd.int_subj = 6, sd_eps = 2) {
#   design <-
#     fixed.factor("obs", levels=c(1:n.obs)) +
#     fixed.factor("cond", levels=c(0:1)) +
#     random.factor("subj", instances=n.subj)
#   dat <- design.codes(design)
#   dat$y <- simLMM(formula = ~ 1 + obs + cond + (1|subj),
#                      data = dat,
#                      Fixef = c(b0, beta_age, beta_cond),
#                      VC_sd = list(c(sd.int_subj), sd_eps),
#                      empirical = TRUE,
#                      verbose = FALSE)
#   return(dat)
# }
# 
# sim_data_1RI.1RS <- function(n.subj = 10, n.obs = 6, n.group = 2, b0 = 10, beta_age = 0, beta_cond = 5, beta_group = 10, sd.int_subj = 6, sd.slope_group = 10, corr_subj.group = .5, sd_eps = 2) {
#   design <-
#     fixed.factor("obs", levels=c(1:n.obs)) +
#     fixed.factor("cond", levels=c(0:1)) +
#     random.factor("subj", instances=n.subj) +
#     random.factor("group", instances = n.group)
#   dat <- design.codes(design)
#   dat$y <- simLMM(formula = ~ 1 + obs + cond + group + (1 + group|subj),
#                      data = dat,
#                      Fixef = c(b0, beta_age, beta_cond, beta_group),
#                      VC_sd = list(c(sd.int_subj, sd.slope_group), sd_eps),
#                      CP = corr_subj.group,
#                      empirical = TRUE,
#                      verbose = FALSE)
#   return(dat)
# }
# 
# sim_data_1RI.2RS <- function(n.subj = 10, n.obs = 6, n.group1 = 2, n.group2 = 2, b0 = 10, beta_age = 0, beta_cond = 5, beta_group1 = 10, beta_group2 = 10, sd.int_subj = 6, sd.slope_group1 = 10, sd.slope_group2 = 15, corr_subj.group = .5, sd_eps = 2) {
#   design <-
#     fixed.factor("obs", levels=c(1:n.obs)) +
#     fixed.factor("cond", levels=c(0:1)) +
#     random.factor("subj", instances=n.subj) +
#     random.factor("group1", instances = n.group1) +
#     random.factor("group2", instances = n.group2)
#   dat <- design.codes(design)
#   dat$y <- simLMM(formula = ~ 1 + obs + cond + group1 + group2 + (1 + group1 + group2|subj),
#                      data = dat,
#                      Fixef = c(b0, beta_age, beta_cond, beta_group1, beta_group2),
#                      VC_sd = list(c(sd.int_subj, sd.slope_group1, sd.slope_group2), sd_eps),
#                      CP = corr_subj.group,
#                      empirical = TRUE,
#                      verbose = FALSE)
#   return(dat)
# }

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

##LRT:
test_lrtstat <- function(data, m.full, m.null, REML = TRUE) {
  full <- lmer(m.full, data = data, REML = REML)
  null <- lmer(m.null, data = data, REML = REML)
  return(pchisq(as.numeric(2 * (logLik(full) - logLik(null))), df = 1, lower = FALSE))
}
m.full <- y ~ age + (1|observation)
m.null <- y ~ (1|observation)
data <- sim_data_1RI(beta_age = .1)
data <- sim_data(model = 1)
REML = FALSE
full <- lmer(m.full, data = data, REML = REML)
null <- lmer(m.null, data = data, REML = REML)
pchisq(as.numeric(2 * (logLik(full) - logLik(null))), df = 1, lower = FALSE)


##t-as-z
test_TasZ.fixed <- function(data, m.full, REML = TRUE) {
  return(summary(lmer(model, data = data, REML = REML))$coefficients[2,4])
}

##KR, SW
#ddf ... Art der Approximation
test_approx.fixed <- function(data, model, REML = TRUE, ddf = "Satterthwaite") {
  return(anova(lmer(model, data = data, REML = REML), ddf = ddf)$`Pr(>F)`[1])
}

##Funktion zur Ausgabe des p-Wertes des fixed effects via parametric bootstrap (package afex)
#mixed auf afex (nutzt pbmodcomp)
#nsim.pb bestimmt anzahl an bootstrap-simulationen von pbmodcomp
#cl erlaubt multicore nutzung (via package parallel)
test_PB.fixed <- function(mode, data, nsim.pb = 1000, cl = NULL) {
  return(suppressMessages(mixed(model, data = data, method = "PB", progress = FALSE, cl = cl, args_test = list(nsim = nsim.pb, cl = cl))$anova_table$`Pr(>PB)`[1]))
}
#suppressMessages: "mixed" will throw a message if numerical variables are not centered on 0, as main effects (of other variables then the numeric one) can be hard to interpret if numerical variables appear in interactions. See Dalal & Zickar (2012).
#kleine Tests haben ergeben, dass es am effizientesten ist, fürs fitten und fürs bootstrap multicore zu nutzen


#Parameter für Simulationen
nsim <- 500
beta_age <- 0 #auf diesen fixed effect wird jeweils getestet
models <- 1:3

#sapply
plan("multisession", workers = detectCores())

#Parameter für parametric bootstrap
nsim.mixed <- 5 #niedriger, weil pro iteration auch noch gebootstrapped wird (mit nsim.pb)
nsim.pb <- 5

###LRT
models_LRT <- c(y ~ 1 + age + (1|observation), y ~ 1 + age + school + (school|observation), y ~ 1 + age + school + class + (school + class|observation))
##REML (nicht empfohlen)
data_LRT.REML_1 <- t(future_replicate(nsim, test_lrtstat(sim_data(model = 1, beta_age = beta_age), models_LRT[[1]], update.formula(models_LRT[[1]],  ~ . - obs))))
colnames(data_LRT.REML_1) <- 1:nsim
data_LRT.REML_2 <- t(future_replicate(nsim, test_lrtstat(sim_data(model = 2, beta_age = beta_age), models_LRT[[2]], update.formula(models_LRT[[2]],  ~ . - obs))))
colnames(data_LRT.REML_2) <- 1:nsim
data_LRT.REML_3 <- t(future_replicate(nsim, test_lrtstat(sim_data(model = 3, beta_age = beta_age), models_LRT[[3]], update.formula(models_LRT[[3]],  ~ . - obs))))
colnames(data_LRT.REML_3) <- 1:nsim
data_LRT.REML_long <- rbind(data_LRT.REML_1, data_LRT.REML_2, data_LRT.REML_3)
data_LRT.REML_long <- as.data.frame(cbind(models, data_LRT.REML_long))
data_LRT.REML_long <- gather(data_LRT.REML_long, sim, p.LRT.REML, 2:ncol(data_LRT.REML_long))

data_LRT.REML_long %>% 
  group_by(model) %>% 
  summarize(prop_LRT.REML = mean(p.LRT.REML <= .05))

##Daten für Plot:
p_LRT.REML <- data_LRT.REML_long %>% 
  group_by(models) %>% 
  summarize(k = sum(p.LRT.REML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(models, p, p_l, p_u) %>% 
  mutate(REML = 1,
         method = 1)

##ML
data_LRT.ML_1 <- t(future_replicate(nsim, test_lrtstat(sim_data(model = 1, beta_age = beta_age), models_LRT[[1]], update.formula(models_LRT[[1]],  ~ . - obs), REML = FALSE)))
colnames(data_LRT.ML_1) <- 1:nsim
data_LRT.ML_2 <- t(future_replicate(nsim, test_lrtstat(sim_data(model = 2, beta_age = beta_age), models_LRT[[2]], update.formula(models_LRT[[2]],  ~ . - obs), REML = FALSE)))
colnames(data_LRT.ML_2) <- 1:nsim
data_LRT.ML_3 <- t(future_replicate(nsim, test_lrtstat(sim_data(model = 3, beta_age = beta_age), models_LRT[[3]], update.formula(models_LRT[[3]],  ~ . - obs), REML = FALSE)))
colnames(data_LRT.ML_3) <- 1:nsim
data_LRT.ML_long <- rbind(data_LRT.ML_1, data_LRT.ML_2, data_LRT.ML_3)
data_LRT.ML_long <- as.data.frame(cbind(model = 1:length(models), data_LRT.ML_long))
data_LRT.ML_long <- gather(data_LRT.ML_long, sim, p.LRT.ML, 2:ncol(data_LRT.ML_long))

data_LRT.ML_long %>% 
  group_by(model) %>% 
  summarize(prop_LRT.ML = mean(p.LRT.ML <= .05))

p_LRT.ML <- data_LRT.ML_long %>% 
  group_by(models) %>% 
  summarize(k = sum(p.LRT.ML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(models, p, p_l, p_u) %>% 
  mutate(ML = 0,
         method = 1)

###t-as-z
##REML
data_TasZ.REML <- t(sapply(models, function(x) future_replicate(nsim, test_TasZ.fixed(sim_data(model = x, beta_age = beta_age), model))))
colnames(data_TasZ.REML) <- 1:nsim
data_TasZ.REML_long <- as.data.frame(cbind(models, data_TasZ.REML))
data_TasZ.REML_long <- gather(data_TasZ.REML_long, sim, p.TasZ.REML, 2:ncol(data_TasZ.REML_long))

data_TasZ.REML_long %>% 
  group_by(models) %>% 
  summarize(prop_TasZ.REML = mean(abs(p.TasZ.REML) >= 1.96))

p_TasZ.REML <- data_TasZ.REML_long %>% 
  group_by(models) %>% 
  summarize(k = sum(abs(p.TasZ.REML) >= 1.96) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(models, p, p_l, p_u) %>% 
  mutate(REML = 1, 
         method = 2)

##ML
data_TasZ.ML <- t(sapply(models, function(x) future_replicate(nsim, test_TasZ.fixed(sim_data_int(beta_age = x), model, REML = FALSE))))
colnames(data_TasZ.ML) <- 1:nsim
data_TasZ.ML_long <- as.data.frame(cbind(models, data_TasZ.ML))
data_TasZ.ML_long <- gather(data_TasZ.ML_long, sim, p.TasZ.ML, 2:ncol(data_TasZ.ML_long))

data_TasZ.ML_long %>% 
  group_by(models) %>% 
  summarize(prop_TasZ.ML = mean(abs(p.TasZ.ML) >= 1.96))

p_TasZ.ML <- data_TasZ.ML_long %>% 
  group_by(models) %>% 
  summarize(k = sum(abs(p.TasZ.ML) >= 1.96) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(models, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 2)

###KR, SW
#anova aus lmertest

##Sattherthwaire, REML
ddf <- "Satterthwaite"
REML <- TRUE
data_SW.REML <- t(sapply(models, function(x) future_replicate(nsim, test_approx.fixed(sim_data_int(beta_age = x), model, REML = REML, ddf = ddf))))
colnames(data_SW.REML) <- 1:nsim
data_SW.REML_long <- as.data.frame(cbind(models, data_SW.REML))
data_SW.REML_long <- gather(data_SW.REML_long, sim, p.SW.REML, 2:ncol(data_SW.REML_long))

data_SW.REML_long %>% 
  group_by(models) %>% 
  summarize(prop_SW.REML = mean(p.SW.REML <= .05))

p_SW.REML <- data_SW.REML_long %>% 
  group_by(models) %>% 
  summarize(k = sum(p.SW.REML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(models, p, p_l, p_u) %>% 
  mutate(REML = 1,
         method = 3)

##Kenward-Roger, REML
ddf <- "Kenward-Roger"
REML <- TRUE
data_KR.REML <- t(sapply(models, function(x) future_replicate(nsim, test_approx.fixed(sim_data_int(beta_age = x), model, REML = TRUE, ddf = "Satterthwaite"))))
colnames(data_KR.REML) <- 1:nsim
data_KR.REML_long <- as.data.frame(cbind(models, data_KR.REML))
data_KR.REML_long <- gather(data_KR.REML_long, sim, p.KR.REML, 2:ncol(data_KR.REML_long))

data_KR.REML_long %>% 
  group_by(models) %>% 
  summarize(prop_KR.REML = mean(p.KR.REML <= .05))

p_KR.REML <- data_KR.REML_long %>% 
  group_by(models) %>% 
  summarize(k = sum(p.KR.REML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(models, p, p_l, p_u) %>% 
  mutate(REML = 1, 
         method = 4)

##Sattherthwaire, ML
ddf <- "Satterthwaite"
REML <- FALSE
data_SW.ML <- t(sapply(models, function(x) future_replicate(nsim, test_approx.fixed(sim_data_int(beta_age = x), model, REML = TRUE, ddf = "Satterthwaite"))))
colnames(data_SW.ML) <- 1:nsim
data_SW.ML_long <- as.data.frame(cbind(models, data_SW.ML))
data_SW.ML_long <- gather(data_SW.ML_long, sim, p.SW.ML, 2:ncol(data_SW.ML_long))

data_SW.ML_long %>% 
  group_by(models) %>% 
  summarize(prop_SW.ML = mean(p.SW.ML <= .05))

p_SW.ML <- data_SW.ML_long %>% 
  group_by(models) %>% 
  summarize(k = sum(p.SW.ML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(models, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 3)

#Kenward-Roger nur für ML möglich!


###parametric bootstrap (nur ML)

#Cluster festlegen (sapply funktioniert nicht)
nc <- detectCores() # number of cores
cl <- makeCluster(rep("localhost", nc)) # make cluster

data_alpha.nB <- t(sapply(models, function(x) replicate(nsim.mixed, test_PB.fixed(model, data = sim_data_int(beta_age = x), nsim.pb = nsim.pb, cl = cl))))
colnames(data_alpha.nB) <- 1:nsim
data_alpha.nB_long <- as.data.frame(cbind(models, data_alpha.nB))
data_alpha.nB_long <- gather(data_alpha.nB_long, sim, p.PB, 2:ncol(data_alpha.nB_long))

data_alpha.nB_long %>% 
  group_by(models) %>% 
  summarize(prop_PB = mean(p.PB <= .05))

p_PB <- data_alpha.nB_long %>% 
  group_by(models) %>% 
  summarize(k = sum(p.PB < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(models, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 5)

### Grafiken der Ergebnisse
data_alpha.n <- rbind(p_TasZ.ML, p_TasZ.REML, p_LRT.ML, p_LRT.REML, p_SW.ML, p_SW.REML, p_KR.REML, p_PB)
data_alpha.n$n.obs <- as.factor(data_alpha.n$n.obs)
data_alpha.n$n.subj <- as.factor(data_alpha.n$n.subj)
data_alpha.n$REML <- factor(data_alpha.n$REML, labels = c("ML", "REML"))
data_alpha.n$method <- factor(data_alpha.n$method, labels = c("LRT", "t-as-z", "Satterthwaite", "Kenward-Roger", "Parametric Bootstrap"))

#alle Methoden
ggplot(data_alpha.n, aes(x = as.factor(n.obs), y = p, col = REML, shape = method)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .05) +
  facet_wrap(~n.subj) +
  ylim(0, .12)

#nur SW und KR
data_alpha.n %>% 
  filter(method %in% c("Satterthwaite", "Kenward-Roger")) %>% 
  ggplot(aes(x = as.factor(n.obs), y = p, col = REML, shape = method)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .05) +
  facet_wrap(~n.subj, nrow = 1) +
  ylim(0, .1)

#nur SW
data_alpha.n %>% 
  filter(method %in% c("Satterthwaite")) %>% 
  ggplot(aes(x = as.factor(n.obs), y = p, col = REML, shape = method)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .05) +
  facet_wrap(~n.subj, nrow = 1) +
  ylim(0, .1)

#nur ML
data_alpha.n %>% 
  filter(REML == "ML") %>% 
  ggplot(aes(x = as.factor(n.obs), y = p, col = method)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .05) +
  facet_wrap(~n.subj, nrow = 1) +
  ylim(0, .12)

#nur REML
data_alpha.n %>% 
  filter(REML == "REML") %>% 
  ggplot(aes(x = as.factor(n.obs), y = p, col = method)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .05) +
  facet_wrap(~n.subj, nrow = 1) +
  ylim(0, .1)

#nur t as z
data_alpha.n %>% 
  filter(method == "t-as-z") %>% 
  ggplot(aes(x = as.factor(n.obs), y = p, col = REML)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .05) +
  facet_wrap(~n.subj, nrow = 1) +
  ylim(0, .12)



library(designr)



summary(lmer(y ~ 1 + obs + cond + group + (1 + group| subj), data = dat))

