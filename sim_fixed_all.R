
path <- ""
n.cores_no <- 0 #number of cores not to be used for parallelization

library(future.apply)
library(tidyverse)
library(lme4)
library(lmerTest)
library(foreach)
library(doRNG)


nsim <- 2
nsim.pb <- 2 #anzahl an simulationen fürs parametric bootstrap
n.bs <- 2 #anzahl an bootstrap-ziehungen

set.seed(1996)

##Datengeneration
#einfaches Modell nur mit random intercept
# y = b0 + b2*cond + (1|subj) + epsilon
sim_data_int <- function(n.subj = 10, n.obs = 6, b0 = 10, beta_cond = 0, sd.int_subj = 6, sd_eps = 2) {
  subj <- rep(1:n.subj, each = n.obs * 2)
  cond <- rep(c(0,1), n.subj*n.obs)
  subj_int <- rep(rnorm(n.subj, 0, sd.int_subj), each = n.obs*2)
  y <- 10 + beta_cond * cond + subj_int + rnorm(length(subj), 0, 2)
  return(data.frame(subj, cond, y))
}

##LRT:
test_lrtstat <- function(data, m.full, m.null, REML = TRUE) {
  full <- lmer(m.full, data = data, REML = REML)
  null <- lmer(m.null, data = data, REML = REML)
  return(pchisq(as.numeric(2 * (logLik(full) - logLik(null))), 1, lower = FALSE))
}

##t-as-z
test_TasZ.fixed <- function(data, m.full, REML = TRUE) {
  return(summary(lmer(model, data = data, REML = REML))$coefficients[2,4])
}

##KR, SW
#ddf ... Art der Approximation
test_approx.fixed <- function(data, model, REML = TRUE, ddf = "Satterthwaite") {
  return(anova(lmer(model, data = data, REML = REML), ddf = ddf)$`Pr(>F)`[1])
}

##Funktion fürs parametric bootstrapping
#benötigt package foreach
#mixed (bz.w pbmodcomp) parallelisiert nicht
#doRNG notwendig, um richtig seed zu setzen
test_PB <- function(data, m.full, m.null, n.bs, REML = FALSE) {
  nullmod <- lmer(m.null, data = data, REML = REML)
  fullmod <- lmer(m.full, data = data, REML = REML)
  lrstat <- numeric(n.bs)
  lrstat <- foreach(i = 1:n.bs, .combine = "c") %dopar% {
    data$y_sim <- unlist(simulate(nullmod))
    null <- lmer(update(m.null, y_sim ~ .), data = data, REML = REML)
    alt <- lmer(update(m.full, y_sim ~ .), data = data, REML = REML)
    as.numeric(2*(logLik(alt)-logLik(null)))
  }
  return(mean(lrstat > as.numeric(2*(logLik(fullmod)-logLik(nullmod)))))
}

#full and null model (LRT):
m.full <- y ~ cond + (1|subj)
m.null <- y ~ (1|subj)

#model
model <- y ~ cond + (1|subj)

#Parameter für Simulationen
beta_cond <- 0 #auf diesen fixed effect wird jeweils getestet
n.subj <- c(4, 6, 10, 16)
n.obs <- c(4, 6, 10, 16)
grid <- expand.grid(n.subj, n.obs)
colnames(grid) <- c("n.subj", "n.obs")

#future_apply
plan("multisession", workers = availableCores() - n.cores_no, gc = TRUE)

###LRT
##REML (nicht empfohlen)
data_LRT.REML <- t(future_apply(grid, 1, function(x) replicate(nsim, test_lrtstat(sim_data_int(n.subj = x[1], n.obs = x[2], beta_cond = beta_cond), m.full, m.null)), future.seed = TRUE))
data_LRT.REML_long <- cbind(grid, data_LRT.REML)
data_LRT.REML_long <- gather(data_LRT.REML_long, sim, p.LRT.REML, (ncol(grid)+1):ncol(data_LRT.REML_long))

p_LRT.REML <- data_LRT.REML_long %>% 
  group_by(n.subj, n.obs) %>% 
  summarize(k = sum(p.LRT.REML < .05) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(n.obs, n.subj, p, p_l, p_u) %>% 
  mutate(REML = 1,
         method = 1)

##ML
data_LRT.ML <- t(future_apply(grid, 1, function(x) replicate(nsim, test_lrtstat(sim_data_int(n.subj = x[1], n.obs = x[2], beta_cond = beta_cond), m.full, m.null, REML = FALSE)), future.seed = TRUE))
data_LRT.ML_long <- cbind(grid, data_LRT.ML)
data_LRT.ML_long <- gather(data_LRT.ML_long, sim, p.LRT.ML, (ncol(grid)+1):ncol(data_LRT.ML_long))

p_LRT.ML <- data_LRT.ML_long %>% 
  group_by(n.subj, n.obs) %>% 
  summarize(k = sum(p.LRT.ML < .05) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(n.obs, n.subj, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 1)

###t-as-z
##REML
data_TasZ.REML <- t(future_apply(grid, 1, function(x) replicate(nsim, test_TasZ.fixed(sim_data_int(n.subj = x[1], n.obs = x[2], beta_cond = beta_cond), model)), future.seed = TRUE))
data_TasZ.REML_long <- cbind(grid, data_TasZ.REML)
data_TasZ.REML_long <- gather(data_TasZ.REML_long, sim, p.TasZ.REML, (ncol(grid)+1):ncol(data_TasZ.REML_long))


p_TasZ.REML <- data_TasZ.REML_long %>% 
  group_by(n.subj, n.obs) %>% 
  summarize(k = sum(abs(p.TasZ.REML) >= qnorm(.975)) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(n.obs, n.subj, p, p_l, p_u) %>% 
  mutate(REML = 1,
         method = 2)

##ML
data_TasZ.ML <- t(future_apply(grid, 1, function(x) replicate(nsim, test_TasZ.fixed(sim_data_int(n.subj = x[1], n.obs = x[2], beta_cond = beta_cond), model, REML = FALSE)), future.seed = TRUE))
data_TasZ.ML_long <- cbind(grid, data_TasZ.ML)
data_TasZ.ML_long <- gather(data_TasZ.ML_long, sim, p.TasZ.ML, (ncol(grid)+1):ncol(data_TasZ.ML_long))

p_TasZ.ML <- data_TasZ.ML_long %>% 
  group_by(n.subj, n.obs) %>% 
  summarize(k = sum(abs(p.TasZ.ML) >= qnorm(.975)) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(n.obs, n.subj, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 2)

###KR, SW
#anova aus lmertest

##Sattherthwaire, REML
ddf <- "Satterthwaite"
REML <- TRUE
data_SW.REML <- t(future_apply(grid, 1, function(x) replicate(nsim, test_approx.fixed(sim_data_int(n.subj = x[1], n.obs = x[2], beta_cond = beta_cond), model, REML = REML, ddf = ddf)), future.seed = TRUE))
data_SW.REML_long <- cbind(grid, data_SW.REML)
data_SW.REML_long <- gather(data_SW.REML_long, sim, p.SW.REML, (ncol(grid)+1):ncol(data_SW.REML_long))

p_SW.REML <- data_SW.REML_long %>% 
  group_by(n.subj, n.obs) %>% 
  summarize(k = sum(p.SW.REML < .05) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(n.obs, n.subj, p, p_l, p_u) %>% 
  mutate(REML = 1,
         method = 3)

##Kenward-Roger, REML
ddf <- "Kenward-Roger"
REML <- TRUE
data_KR.REML <- t(future_apply(grid, 1, function(x) replicate(nsim, test_approx.fixed(sim_data_int(n.subj = x[1], n.obs = x[2], beta_cond = beta_cond), model, REML = REML, ddf = ddf)), future.seed = TRUE))
data_KR.REML_long <- cbind(grid, data_KR.REML)
data_KR.REML_long <- gather(data_KR.REML_long, sim, p.KR.REML, (ncol(grid)+1):ncol(data_KR.REML_long))

p_KR.REML <- data_KR.REML_long %>% 
  group_by(n.subj, n.obs) %>% 
  summarize(k = sum(p.KR.REML < .05) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(n.obs, n.subj, p, p_l, p_u) %>% 
  mutate(REML = 1,
         method = 4)

##Sattherthwaire, ML
ddf <- "Satterthwaite"
REML <- FALSE
data_SW.ML <- t(future_apply(grid, 1, function(x) replicate(nsim, test_approx.fixed(sim_data_int(n.subj = x[1], n.obs = x[2], beta_cond = beta_cond), model, REML = REML, ddf = ddf)), future.seed = TRUE))
data_SW.ML_long <- cbind(grid, data_SW.ML)
data_SW.ML_long <- gather(data_SW.ML_long, sim, p.SW.ML, (ncol(grid)+1):ncol(data_SW.ML_long))

p_SW.ML <- data_SW.ML_long %>% 
  group_by(n.subj, n.obs) %>% 
  summarize(k = sum(p.SW.ML < .05) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(n.obs, n.subj, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 3)

#Kenward-Roger nur für ML möglich!

#future clusters stoppen
future:::ClusterRegistry("stop")


###parametric bootstrap
##ML
REML = FALSE
#Cluster festlegen
n.cores <- parallel::detectCores() - n.cores_no
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)

#set seed für foreach
registerDoRNG(1996)

data_PB.ML <- t(apply(grid, 1, function(x) replicate(nsim.pb, test_PB(data = sim_data_int(n.subj = x[1], n.obs = x[2], beta_cond = beta_cond), m.full = m.full,m.null = m.null, n.bs = n.bs, REML = REML))))
data_PB.ML_long <- cbind(grid, data_PB.ML)
data_PB.ML_long <- gather(data_PB.ML_long, sim, p.PB.ML, (ncol(grid)+1):ncol(data_PB.ML_long))

p_PB.ML <- data_PB.ML_long %>% 
  group_by(n.subj, n.obs) %>% 
  summarize(k = sum(p.PB.ML < .05) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(n.obs, n.subj, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 5)

##REML
REML = TRUE

#set seed für foreach
registerDoRNG(123)

data_PB.REML <- t(apply(grid, 1, function(x) replicate(nsim.pb, test_PB(data = sim_data_int(n.subj = x[1], n.obs = x[2], beta_cond = beta_cond), m.full = m.full, m.null = m.null ,n.bs = n.bs, REML = REML))))
parallel::stopCluster(cl = my.cluster)
data_PB.REML_long <- cbind(grid, data_PB.REML)
data_PB.REML_long <- gather(data_PB.REML_long, sim, p.PB.REML, (ncol(grid)+1):ncol(data_PB.REML_long))

p_PB.REML <- data_PB.REML_long %>% 
  group_by(n.subj, n.obs) %>% 
  summarize(k = sum(p.PB.REML < .05) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(n.obs, n.subj, p, p_l, p_u) %>% 
  mutate(REML = 1,
         method = 5)

#gesamte Daten
data_n <- rbind(p_TasZ.ML, p_TasZ.REML, p_LRT.ML, p_LRT.REML, p_SW.ML, p_SW.REML, p_KR.REML, p_PB.ML, p_PB.REML)
data_n$n.obs <- as.factor(data_n$n.obs)
data_n$n.subj <- as.factor(data_n$n.subj)
data_n$REML <- factor(data_n$REML, labels = c("ML", "REML"))
data_n$method <- factor(data_n$method, labels = c("LRT", "t-as-z", "Satterthwaite", "Kenward-Roger", "Parametric Bootstrap"))


save(data_n, file = paste0(path, "data_n.RData"))

#####################################################################################################################################################################################################
#Missing Values
#####################################################################################################################################################################################################

##missing values erzeugen (MCAR)
missing.val <- function(data, p) {
  return(data[-sample(1:nrow(data), size = p * nrow(data), replace = FALSE),])
}

#full and null model (LRT):
m.full <- y ~ cond + (1|subj)
m.null <- y ~ (1|subj)

#model
model <- y ~ cond + (1|subj)

#Parameter für Simulationen
beta_cond <- 0 #auf diesen fixed effect wird jeweils getestet
p.missing <- c(.1, .3, .5)

set.seed(1996)

plan("multisession", workers = detectCores() - n.cores_no)

###LRT
##REML (nicht empfohlen)
data_LRT.REML <- t(sapply(p.missing, function(x) future_replicate(nsim, test_lrtstat(missing.val(sim_data_int(n.subj = 100), p = x), m.full, m.null))))
colnames(data_LRT.REML) <- 1:nsim
data_LRT.REML_long <- as.data.frame(cbind(p.missing, data_LRT.REML))
data_LRT.REML_long <- gather(data_LRT.REML_long, sim, p.LRT.REML, 2:ncol(data_LRT.REML_long))

##Daten für Plot:
p_LRT.REML <- data_LRT.REML_long %>% 
  group_by(p.missing) %>% 
  summarize(k = sum(p.LRT.REML < .05) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(p.missing, p, p_l, p_u) %>% 
  mutate(REML = 1,
         method = 1)

##ML
data_LRT.ML <- t(sapply(p.missing, function(x) future_replicate(nsim, test_lrtstat(missing.val(sim_data_int(), p = x), m.full, m.null, REML = FALSE))))
colnames(data_LRT.ML) <- 1:nsim
data_LRT.ML_long <- as.data.frame(cbind(p.missing, data_LRT.ML))
data_LRT.ML_long <- gather(data_LRT.ML_long, sim, p.LRT.ML, 2:ncol(data_LRT.ML_long))

p_LRT.ML <- data_LRT.ML_long %>% 
  group_by(p.missing) %>% 
  summarize(k = sum(p.LRT.ML < .05) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(p.missing, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 1)

###t-as-z
##REML
data_TasZ.REML <- t(sapply(p.missing, function(x) future_replicate(nsim, test_TasZ.fixed(missing.val(sim_data_int(), p = x), model))))
colnames(data_TasZ.REML) <- 1:nsim
data_TasZ.REML_long <- as.data.frame(cbind(p.missing, data_TasZ.REML))
data_TasZ.REML_long <- gather(data_TasZ.REML_long, sim, p.TasZ.REML, 2:ncol(data_TasZ.REML_long))

p_TasZ.REML <- data_TasZ.REML_long %>% 
  group_by(p.missing) %>% 
  summarize(k = sum(abs(p.TasZ.REML) >= qnorm(.975)) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(p.missing, p, p_l, p_u) %>% 
  mutate(REML = 1, 
         method = 2)

##ML
data_TasZ.ML <- t(sapply(p.missing, function(x) future_replicate(nsim, test_TasZ.fixed(missing.val(sim_data_int(), p = x), model, REML = FALSE))))
colnames(data_TasZ.ML) <- 1:nsim
data_TasZ.ML_long <- as.data.frame(cbind(p.missing, data_TasZ.ML))
data_TasZ.ML_long <- gather(data_TasZ.ML_long, sim, p.TasZ.ML, 2:ncol(data_TasZ.ML_long))

p_TasZ.ML <- data_TasZ.ML_long %>% 
  group_by(p.missing) %>% 
  summarize(k = sum(abs(p.TasZ.ML) >= qnorm(.975)) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(p.missing, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 2)

###KR, SW
#anova aus lmertest

##Sattherthwaire, REML
ddf <- "Satterthwaite"
REML <- TRUE
data_SW.REML <- t(sapply(p.missing, function(x) future_replicate(nsim, test_approx.fixed(missing.val(sim_data_int(), p = x), model, REML = REML, ddf = ddf))))
colnames(data_SW.REML) <- 1:nsim
data_SW.REML_long <- as.data.frame(cbind(p.missing, data_SW.REML))
data_SW.REML_long <- gather(data_SW.REML_long, sim, p.SW.REML, 2:ncol(data_SW.REML_long))

p_SW.REML <- data_SW.REML_long %>% 
  group_by(p.missing) %>% 
  summarize(k = sum(p.SW.REML < .05) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(p.missing, p, p_l, p_u) %>% 
  mutate(REML = 1,
         method = 3)

##Kenward-Roger, REML
ddf <- "Kenward-Roger"
REML <- TRUE
data_KR.REML <- t(sapply(p.missing, function(x) future_replicate(nsim, test_approx.fixed(missing.val(sim_data_int(), p = x), model, REML = REML, ddf = ddf))))
colnames(data_KR.REML) <- 1:nsim
data_KR.REML_long <- as.data.frame(cbind(p.missing, data_KR.REML))
data_KR.REML_long <- gather(data_KR.REML_long, sim, p.KR.REML, 2:ncol(data_KR.REML_long))

p_KR.REML <- data_KR.REML_long %>% 
  group_by(p.missing) %>% 
  summarize(k = sum(p.KR.REML < .05) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(p.missing, p, p_l, p_u) %>% 
  mutate(REML = 1, 
         method = 4)

##Sattherthwaire, ML
ddf <- "Satterthwaite"
REML <- FALSE
data_SW.ML <- t(sapply(p.missing, function(x) future_replicate(nsim, test_approx.fixed(missing.val(sim_data_int(), p = x), model, REML = REML, ddf = ddf))))
colnames(data_SW.ML) <- 1:nsim
data_SW.ML_long <- as.data.frame(cbind(p.missing, data_SW.ML))
data_SW.ML_long <- gather(data_SW.ML_long, sim, p.SW.ML, 2:ncol(data_SW.ML_long))

p_SW.ML <- data_SW.ML_long %>% 
  group_by(p.missing) %>% 
  summarize(k = sum(p.SW.ML < .05) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(p.missing, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 3)

#Kenward-Roger nur für ML möglich!

#future cluster stoppen
future:::ClusterRegistry("stop")

###parametric bootstrap
##ML
REML = FALSE
#Cluster festlegen
n.cores <- parallel::detectCores() - n.cores_no
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)

#set seed für foreach
registerDoRNG(1996)

data_PB.ML <- t(sapply(p.missing, function(x) replicate(nsim.pb, test_PB(data = missing.val(sim_data_int(), p = x), m.full = m.full, m.null = m.null, n.bs = n.bs, REML = REML))))
colnames(data_PB.ML) <- 1:nsim.pb
data_PB.ML_long <- as.data.frame(cbind(p.missing, data_PB.ML))
data_PB.ML_long <- gather(data_PB.ML_long, sim, p.PB.ML, 2:ncol(data_PB.ML_long))

p_PB <- data_PB.ML_long %>% 
  group_by(p.missing) %>% 
  summarize(k = sum(p.PB.ML < .05) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(p.missing, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 5)

##REML
REML = TRUE

#set seed für foreach
registerDoRNG(1996)

data_PB.REML <- t(sapply(p.missing, function(x) replicate(nsim.pb, test_PB(data = missing.val(sim_data_int(), p = x), m.full = m.full, m.null = m.null, n.bs = n.bs, REML = REML))))
parallel::stopCluster(cl = my.cluster)
colnames(data_PB.REML) <- 1:nsim.pb
data_PB.REML_long <- as.data.frame(cbind(p.missing, data_PB.REML))
data_PB.REML_long <- gather(data_PB.REML_long, sim, p.PB.REML, 2:ncol(data_PB.REML_long))

p_PB <- data_PB.REML_long %>% 
  group_by(p.missing) %>% 
  summarize(k = sum(p.PB.REML < .05) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(p.missing, p, p_l, p_u) %>% 
  mutate(REML = 1,
         method = 5)


data_missing <- rbind(p_TasZ.ML, p_TasZ.REML, p_LRT.ML, p_LRT.REML, p_SW.ML, p_SW.REML, p_KR.REML, p_PB)
data_missing$p.missing <- as.factor(data_missing$p.missing)
data_missing$REML <- factor(data_missing$REML, labels = c("ML", "REML"))
data_missing$method <- factor(data_missing$method, labels = c("LRT", "t-as-z", "Satterthwaite", "Kenward-Roger", "Parametric Bootstrap"))

save(data_missing, file = paste0(path, "data_missing.RData"))

#####################################################################################################################################################################################################
#Unbalanced Design
#####################################################################################################################################################################################################


##Datengeneration
#einfaches Modell nur mit random intercept
# y = b0 + beta_group * group + (1|cond)
#n.subj und n.obs müssen gerade sein
sim_data_int_unb <- function(n.group1 = 10, n.group2 = 10, n.cond = 4, b0 = 10, beta_group = 0, sd.int_cond = 6, sd_eps = 1) {
  group <- c(rep(0, n.group1 * n.cond), rep(1, n.group2 * n.cond))
  cond <- rep(1:n.cond, length(group)/n.cond)
  cond_int <- rep(rnorm(n.cond, 0, sd.int_cond), length(group)/4)
  y <- b0 + beta_group * group + cond_int + rnorm(length(group), 0, sd_eps)
  return(data.frame(group = as.factor(group), cond = as.factor(cond), y))
}

#full and null model (LRT):
m.full <- y ~ group + (1|cond)
m.null <- y ~ (1|cond)

#model
model <- y ~ group + (1|cond)

#Parameter für Simulationen
beta_group <- 0 #auf diesen fixed effect wird jeweils getestet
n.group1 <- c(25, 40, 50, 75, 90, 150) 
n.group2 <- c(25, 10, 50, 25, 10, 50)
grid <- cbind(n.group1, n.group2)
#colnames(grid) <- c("n.group1", "n.group2")

set.seed(1996)

plan("multisession", workers = detectCores() - n.cores_no)

###LRT
##REML (nicht empfohlen)
data_LRT.REML <- t(apply(grid, 1, function(x) future_replicate(nsim, test_lrtstat(data = sim_data_int_unb(n.group1 = x[1], n.group2 = x[2]), m.full = m.full, m.null = m.null))))
colnames(data_LRT.REML) <- 1:nsim
data_LRT.REML_long <- as.data.frame(cbind(grid, data_LRT.REML))
data_LRT.REML_long <- gather(data_LRT.REML_long, sim, p.LRT.REML, (ncol(grid)+1):ncol(data_LRT.REML_long))

##Daten für Plot:
p_LRT.REML <- data_LRT.REML_long %>% 
  group_by(n.group1, n.group2) %>% 
  summarize(k = sum(p.LRT.REML < .05) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(n.group1, n.group2, p, p_l, p_u) %>% 
  mutate(REML = 1,
         method = 1)

##ML
data_LRT.ML <- t(apply(grid, 1, function(x) future_replicate(nsim, test_lrtstat(sim_data_int_unb(n.group1 = x[1], n.group2 = x[2]), m.full, m.null, REML = FALSE))))
colnames(data_LRT.ML) <- 1:nsim
data_LRT.ML_long <- as.data.frame(cbind(grid, data_LRT.ML))
data_LRT.ML_long <- gather(data_LRT.ML_long, sim, p.LRT.ML, (ncol(grid)+1):ncol(data_LRT.ML_long))

p_LRT.ML <- data_LRT.ML_long %>% 
  group_by(n.group1, n.group2) %>% 
  summarize(k = sum(p.LRT.ML < .05) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(n.group1, n.group2, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 1)

###t-as-z
##REML
data_TasZ.REML <- t(apply(grid, 1, function(x) future_replicate(nsim, test_TasZ.fixed(sim_data_int_unb(n.group1 = x[1], n.group2 = x[2]), model))))
colnames(data_TasZ.REML) <- 1:nsim
data_TasZ.REML_long <- as.data.frame(cbind(grid, data_TasZ.REML))
data_TasZ.REML_long <- gather(data_TasZ.REML_long, sim, p.TasZ.REML, (ncol(grid)+1):ncol(data_TasZ.REML_long))

p_TasZ.REML <- data_TasZ.REML_long %>% 
  group_by(n.group1, n.group2) %>% 
  summarize(k = sum(abs(p.TasZ.REML) >= qnorm(.975)) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(n.group1, n.group2, p, p_l, p_u) %>% 
  mutate(REML = 1, 
         method = 2)

##ML
data_TasZ.ML <- t(apply(grid, 1, function(x) future_replicate(nsim, test_TasZ.fixed(sim_data_int_unb(n.group1 = x[1], n.group2 = x[2]), model, REML = FALSE))))
colnames(data_TasZ.ML) <- 1:nsim
data_TasZ.ML_long <- as.data.frame(cbind(grid, data_TasZ.ML))
data_TasZ.ML_long <- gather(data_TasZ.ML_long, sim, p.TasZ.ML, (ncol(grid)+1):ncol(data_TasZ.ML_long))

p_TasZ.ML <- data_TasZ.ML_long %>% 
  group_by(n.group1, n.group2) %>% 
  summarize(k = sum(abs(p.TasZ.ML) >= qnorm(.975)) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(n.group1, n.group2, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 2)

###KR, SW
#anova aus lmertest

##Sattherthwaire, REML
ddf <- "Satterthwaite"
REML <- TRUE
data_SW.REML <- t(apply(grid, 1, function(x) future_replicate(nsim, test_approx.fixed(sim_data_int_unb(n.group1 = x[1], n.group2 = x[2]), model, REML = REML, ddf = ddf))))
colnames(data_SW.REML) <- 1:nsim
data_SW.REML_long <- as.data.frame(cbind(grid, data_SW.REML))
data_SW.REML_long <- gather(data_SW.REML_long, sim, p.SW.REML, (ncol(grid)+1):ncol(data_SW.REML_long))

p_SW.REML <- data_SW.REML_long %>% 
  group_by(n.group1, n.group2) %>% 
  summarize(k = sum(p.SW.REML < .05) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(n.group1, n.group2, p, p_l, p_u) %>% 
  mutate(REML = 1,
         method = 3)

##Kenward-Roger, REML
ddf <- "Kenward-Roger"
REML <- TRUE
data_KR.REML <- t(apply(grid, 1, function(x) future_replicate(nsim, test_approx.fixed(sim_data_int_unb(n.group1 = x[1], n.group2 = x[2]), model, REML = REML, ddf = ddf))))
colnames(data_KR.REML) <- 1:nsim
data_KR.REML_long <- as.data.frame(cbind(grid, data_KR.REML))
data_KR.REML_long <- gather(data_KR.REML_long, sim, p.KR.REML, (ncol(grid)+1):ncol(data_KR.REML_long))

p_KR.REML <- data_KR.REML_long %>% 
  group_by(n.group1, n.group2) %>% 
  summarize(k = sum(p.KR.REML < .05) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(n.group1, n.group2, p, p_l, p_u) %>% 
  mutate(REML = 1, 
         method = 4)

##Sattherthwaire, ML
ddf <- "Satterthwaite"
REML <- FALSE
data_SW.ML <- t(apply(grid, 1, function(x) future_replicate(nsim, test_approx.fixed(sim_data_int_unb(n.group1 = x[1], n.group2 = x[2]), model, REML = REML, ddf = ddf))))
colnames(data_SW.ML) <- 1:nsim
data_SW.ML_long <- as.data.frame(cbind(grid, data_SW.ML))
data_SW.ML_long <- gather(data_SW.ML_long, sim, p.SW.ML, (ncol(grid)+1):ncol(data_SW.ML_long))

p_SW.ML <- data_SW.ML_long %>% 
  group_by(n.group1, n.group2) %>% 
  summarize(k = sum(p.SW.ML < .05) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(n.group1, n.group2, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 3)

#Kenward-Roger nur für ML möglich!

#future cluster stoppen
future:::ClusterRegistry("stop")

###parametric bootstrap
##ML
REML = FALSE
#Cluster festlegen
n.cores <- parallel::detectCores() - n.cores_no
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)

#set seed für foreach
registerDoRNG(1996)

data_PB.ML <- t(apply(grid, 1, function(x) replicate(nsim.pb, test_PB(data = sim_data_int_unb(n.group1 = x[1], n.group2 = x[2]), m.full = m.full,m.null = m.null, n.bs = n.bs, REML = REML))))
colnames(data_PB.ML) <- 1:nsim.pb
data_PB.ML_long <- as.data.frame(cbind(grid, data_PB.ML))
data_PB.ML_long <- gather(data_PB.ML_long, sim, p.PB.ML, (ncol(grid)+1):ncol(data_PB.ML_long))

p_PB.ML <- data_PB.ML_long %>% 
  group_by(n.group1, n.group2) %>% 
  summarize(k = sum(p.PB.ML < .05) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(n.group1, n.group2, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 5)

##REML
REML = TRUE

#set seed für foreach
registerDoRNG(1996)

data_PB.REML <- t(apply(grid, 1, function(x) replicate(nsim.pb, test_PB(data = sim_data_int_unb(n.group1 = x[1], n.group2 = x[2]), m.full = m.full,m.null = m.null, n.bs = n.bs, REML = REML))))
parallel::stopCluster(cl = my.cluster)
colnames(data_PB.REML) <- 1:nsim.pb
data_PB.REML_long <- as.data.frame(cbind(grid, data_PB.REML))
data_PB.REML_long <- gather(data_PB.REML_long, sim, p.PB.REML, (ncol(grid)+1):ncol(data_PB.REML_long))

p_PB.REML <- data_PB.REML_long %>% 
  group_by(n.group1, n.group2) %>% 
  summarize(k = sum(p.PB.REML < .05) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(n.group1, n.group2, p, p_l, p_u) %>% 
  mutate(REML = 1,
         method = 5)

data_unbalanced <- rbind(p_TasZ.ML, p_TasZ.REML, p_LRT.ML, p_LRT.REML, p_SW.ML, p_SW.REML, p_KR.REML, p_PB.ML, p_PB.REML)
data_unbalanced$n.group1 <- as.factor(data_unbalanced$n.group1)
data_unbalanced$n.group2 <- as.factor(data_unbalanced$n.group2)
data_unbalanced$REML <- factor(data_unbalanced$REML, labels = c("ML", "REML"))
data_unbalanced$method <- factor(data_unbalanced$method, labels = c("LRT", "t-as-z", "Satterthwaite", "Kenward-Roger", "Parametric Bootstrap"))

save(data_unbalanced, file = paste0(path, "data_unbalanced.RData"))

#####################################################################################################################################################################################################
#Effect Size
#####################################################################################################################################################################################################

#full and null model (LRT):
m.full <- y ~ cond + (1|subj)
m.null <- y ~ (1|subj)

#model
model <- y ~ cond + (1|subj)

#Parameter für Simulationen
beta_cond <- 0 #auf diesen fixed effect wird jeweils getestet
ES <- seq(0, 1.4, .2)
n.obs <- c(4, 6, 10, 16)
grid <- expand.grid(ES, n.obs)
colnames(grid) <- c("ES", "n.obs")

set.seed(1996)

#future_replicate
plan("multisession", workers = detectCores() - n.cores_no)

###LRT
##REML (nicht empfohlen)
data_LRT.REML <- t(apply(grid, 1, function(x) replicate(nsim, test_lrtstat(sim_data_int(beta_cond = x[1], n.obs = x[2]), m.full, m.null))))
colnames(data_LRT.REML) <- 1:nsim
data_LRT.REML_long <- as.data.frame(cbind(grid, data_LRT.REML))
data_LRT.REML_long <- gather(data_LRT.REML_long, sim, p.LRT.REML, (ncol(grid)+1):ncol(data_LRT.REML_long))

##Daten für Plot:
p_LRT.REML <- data_LRT.REML_long %>% 
  group_by(ES, n.obs) %>% 
  summarize(k = sum(p.LRT.REML < .05) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(ES, n.obs, p, p_l, p_u) %>% 
  mutate(REML = 1,
         method = 1)

##ML
data_LRT.ML <- t(apply(grid, 1, function(x) replicate(nsim, test_lrtstat(sim_data_int(beta_cond = x[1], n.obs = x[2]), m.full, m.null, REML = FALSE))))
data_LRT.ML_long <- as.data.frame(cbind(grid, data_LRT.ML))
data_LRT.ML_long <- gather(data_LRT.ML_long, sim, p.LRT.ML, (ncol(grid)+1):ncol(data_LRT.ML_long))

p_LRT.ML <- data_LRT.ML_long %>% 
  group_by(ES, n.obs) %>% 
  summarize(k = sum(p.LRT.ML < .05) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(ES, n.obs, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 1)

###t-as-z
##REML
data_TasZ.REML <- t(apply(grid, 1, function(x) replicate(nsim, test_TasZ.fixed(sim_data_int(beta_cond = x[1], x[2]), model))))
data_TasZ.REML_long <- as.data.frame(cbind(grid, data_TasZ.REML))
data_TasZ.REML_long <- gather(data_TasZ.REML_long, sim, p.TasZ.REML, (ncol(grid)+1):ncol(data_TasZ.REML_long))

p_TasZ.REML <- data_TasZ.REML_long %>% 
  group_by(ES, n.obs) %>% 
  summarize(k = sum(abs(p.TasZ.REML) >= qnorm(.975)) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(ES, n.obs, p, p_l, p_u) %>% 
  mutate(REML = 1, 
         method = 2)

##ML
data_TasZ.ML <- t(apply(grid, 1, function(x) replicate(nsim, test_TasZ.fixed(sim_data_int(beta_cond = x[1], n.obs = x[2]), model, REML = FALSE))))
data_TasZ.ML_long <- as.data.frame(cbind(grid, data_TasZ.ML))
data_TasZ.ML_long <- gather(data_TasZ.ML_long, sim, p.TasZ.ML, (ncol(grid)+1):ncol(data_TasZ.ML_long))

p_TasZ.ML <- data_TasZ.ML_long %>% 
  group_by(ES, n.obs) %>% 
  summarize(k = sum(abs(p.TasZ.ML) >= qnorm(.975)) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(ES, n.obs, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 2)

###KR, SW
#anova aus lmertest

##Sattherthwaire, REML
ddf <- "Satterthwaite"
REML <- TRUE
data_SW.REML <- t(apply(grid, 1, function(x) replicate(nsim, test_approx.fixed(sim_data_int(beta_cond = x[1], n.obs = x[2]), model, REML = REML, ddf = ddf))))
data_SW.REML_long <- as.data.frame(cbind(grid, data_SW.REML))
data_SW.REML_long <- gather(data_SW.REML_long, sim, p.SW.REML, (ncol(grid)+1):ncol(data_SW.REML_long))

p_SW.REML <- data_SW.REML_long %>% 
  group_by(ES, n.obs) %>% 
  summarize(k = sum(p.SW.REML < .05) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(ES, n.obs, p, p_l, p_u) %>% 
  mutate(REML = 1,
         method = 3)

##Kenward-Roger, REML
ddf <- "Kenward-Roger"
REML <- TRUE
data_KR.REML <- t(apply(grid, 1, function(x) replicate(nsim, test_approx.fixed(sim_data_int(beta_cond = x[1], n.obs = x[2]), model, REML = REML, ddf = ddf))))
data_KR.REML_long <- as.data.frame(cbind(grid, data_KR.REML))
data_KR.REML_long <- gather(data_KR.REML_long, sim, p.KR.REML, (ncol(grid)+1):ncol(data_KR.REML_long))

p_KR.REML <- data_KR.REML_long %>% 
  group_by(ES, n.obs) %>% 
  summarize(k = sum(p.KR.REML < .05) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(ES, n.obs, p, p_l, p_u) %>% 
  mutate(REML = 1, 
         method = 4)

##Sattherthwaire, ML
ddf <- "Satterthwaite"
REML <- FALSE
data_SW.ML <- t(apply(grid, 1, function(x) replicate(nsim, test_approx.fixed(sim_data_int(beta_cond = x[1], n.obs = x[2]), model, REML = REML, ddf = ddf))))
data_SW.ML_long <- as.data.frame(cbind(grid, data_SW.ML))
data_SW.ML_long <- gather(data_SW.ML_long, sim, p.SW.ML, (ncol(grid)+1):ncol(data_SW.ML_long))

p_SW.ML <- data_SW.ML_long %>% 
  group_by(ES, n.obs) %>% 
  summarize(k = sum(p.SW.ML < .05) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(ES, n.obs, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 3)

#Kenward-Roger nur für ML möglich!

#future cluster stoppen
future:::ClusterRegistry("stop")

###parametric bootstrap
##ML
REML = FALSE
#Cluster festlegen
n.cores <- parallel::detectCores() - n.cores_no
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)

#set seed für foreach
registerDoRNG(1996)

data_PB.ML <- t(apply(grid, 1, function(x) replicate(nsim.pb, test_PB(data = sim_data_int(beta_cond = x[1], n.obs = x[2]), m.full = m.full, m.null = m.null, n.bs = n.bs, REML = REML))))
data_PB.ML_long <- as.data.frame(cbind(grid, data_PB.ML))
data_PB.ML_long <- gather(data_PB.ML_long, sim, p.PB.ML, (ncol(grid)+1):ncol(data_PB.ML_long))

p_PB <- data_PB.ML_long %>% 
  group_by(ES) %>% 
  summarize(k = sum(p.PB.ML < .05) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(ES, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 5)

##REML
REML = TRUE

#set seed für foreach
registerDoRNG(1996)

data_PB.REML <- t(apply(grid, 1, function(x) replicate(nsim.pb, test_PB(data = sim_data_int(beta_cond = x[1], n.obs = x[2]), m.full = m.full, m.null = m.null, n.bs = n.bs, REML = REML))))
parallel::stopCluster(cl = my.cluster)
data_PB.REML_long <- as.data.frame(cbind(grid, data_PB.REML))
data_PB.REML_long <- gather(data_PB.REML_long, sim, p.PB.REML, (ncol(grid)+1):ncol(data_PB.REML_long))

p_PB <- data_PB.REML_long %>% 
  group_by(ES) %>% 
  summarize(k = sum(p.PB.REML < .05) + qnorm(.975)^2/2,
            n = n() + qnorm(.975)^2,
            p = k/n,
            p_l = p - qnorm(.975) * sqrt(p*(1-p)/n),
            p_u = p + qnorm(.975) * sqrt(p*(1-p)/n)) %>% 
  select(ES, p, p_l, p_u) %>% 
  mutate(REML = 1,
         method = 5)

### Grafiken der Ergebnisse
data_ES <- rbind(p_TasZ.ML, p_TasZ.REML, p_LRT.ML, p_LRT.REML, p_SW.ML, p_SW.REML, p_KR.REML, p_PB.ML, p_PB.REML)
data_ES$ES <- as.factor(data_ES$ES)
data_ES$n.obs <- as.factor(data_ES$n.obs)
data_ES$REML <- factor(data_ES$REML, labels = c("ML", "REML"))
data_ES$method <- factor(data_ES$method, labels = c("LRT", "t-as-z", "Satterthwaite", "Kenward-Roger", "Parametric Bootstrap"))

save(data_ES, file = paste0(path, "data_ES.RData"))
