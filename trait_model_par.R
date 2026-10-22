# R version
# R 4.5.2

# Attached packages
# - foreach           (v1.5.2)
# - doParallel        (v1.0.17)


library(foreach)
library(doParallel)

# Detect and define the number of cores to use
cores <- detectCores()
# Register the parallel backend
cl <- makeCluster(cores-4)
registerDoParallel(cl)


# Mutualistic model function
lotka<-function(t,y,parameters){
  with(as.list(c(y,parameters)),{
    Xx<-y
    XP<-Xx[1:M]
    XA<-Xx[(M+1):(M+N)]
    bet<-matrix(Xx[(M+N+1):(M*N+M+N)],M,N)
    xP<-Xx[(M*N+M+N+1):(M*N+2*M+N)]
    yA<-Xx[(M*N+2*M+N+1):(M*N+2*M+2*N)]


    # carrying capacity
    KP=kP*exp(-((xP-xm)^2)/(2*(sd_kP)^2))
    KA=kA*exp(-((yA-ym)^2)/(2*(sd_kA)^2))
    # within guild dissimilarity matrix
    disMatP=matrix(xP,M,M)-matrix(xP,M,M,byrow = T)
    disMatA=matrix(yA,N,N, byrow = T)-matrix(yA,N,N)
    # cross guild dissimilarity matrix
    disMat=c(matrix(xP,M,1))-matrix(yA,M,N, byrow = T)
    # competition kernel
    cP=1*exp(-(disMatP^2)/(2*(sd_cP)^2))
    cA=1*exp(-(disMatA^2)/(2*(sd_cA)^2))

    # benefit function
    sigma_P=c*exp(-(disMat^2)/(2*(sd_m)^2))
    sigma_A=c*exp(-(disMat^2)/(2*(sd_m)^2))

    dXP<-(rP-(rP*(cP%*%XP))/KP+(((sigma_P*bet)%*%XA)*(alpha))/(w+((sigma_A*bet)%*%XA)))

    dXA<-(rA-(rA*(cA%*%XA))/KA+(t(sigma_A*bet)%*%((alpha*XP)/(w+((sigma_A*bet)%*%XA)))))

    dX<-(Xx[1:(M+N)]*rbind(dXP,dXA))

    dbet<-((bet%*%diag(G,ncol = N))*sweep((sigma_A*c(alpha)*c(XP))/c(w+((sigma_A*bet)%*%XA)),2,
                                          colSums((bet*sigma_A*c(alpha)*c(XP))/c(w+((sigma_A*bet)%*%XA))),"-"))


    #  Selection gradient of plant competitive trait
    gx=c(rP*XP/KP)*rowSums(cP*(-(xP-xm)/sd_kP^2+disMatP/sd_cP^2))+(((c(alpha)*bet^2*disMat*sigma_P*sigma_A)/c(sd_m^2*(((bet*sigma_A)%*%XA)+w)^2))%*%XA^2)-(((c(alpha)*bet*disMat*sigma_P)/c(sd_m^2*(((bet*sigma_A)%*%XA)+w)))%*%XA^2)


    # Selection gradient of animal competitive trait
    gy=colSums((cA*(matrix((-(yA-ym)/sd_kA^2),N,N,byrow = TRUE)+disMatA/sd_cA^2))%*%diag(c((rA*XA/KA)),ncol = N))+colSums((c(XP*alpha)*bet*disMat*sigma_A)/c(sd_m^2*(((bet*sigma_A)%*%XA)+w)))

    # canonical equation of adaptive dynamics (trait dynamics)
    dxP=m*XP*gx
    dyA=m*XA*gy

    return(list(c(dX,dbet,dxP,dyA)))
  })
}


rao_quadratic_entropy <- function(n, x){
  # Normalize n to sum to 1 (recommended for diversity indices)
  P <- n / sum(n)

  # Compute pairwise trait distances
  D <- abs(outer(x, x, "-"))

  # Compute Rao's Q
  Q <- sum(outer(P, P) * D)
  return(Q)
}

pielou_evenness <- function(n) {

  # Relative abundances
  p <- n / sum(n)

  # Shannon diversity
  H <- -sum(p * log(p))

  # Maximum possible Shannon for S species
  S <- length(p)
  Hmax <- log(S)

  # Pielou's evenness
  J <- H / Hmax
  return(J)
}

lotka_deriv <- function(y, parameters){
  unlist(lotka(t = NULL, y = y, parameters = parameters))
}

# --- Compute J at final state ---
compute_resilience <- function(y_final, parameters){
  J <- jacobian(func = function(y) lotka_deriv(y, parameters),
                x    = y_final)
  eig <- eigen(J)
  lambda1 <- max(Re(eig$values))
  resilience <- -lambda1
  return(resilience)

}

run_scenario <- function(M, N, lotka, X0, bet0, xP0, yA0,
                         rP, rA, alpha, w, c,
                         kP, kA, xm, ym,
                         sd_m, sd_kP, sd_kA, sd_cP, sd_cA,
                         G, m, suffix = "", timeout = 30) {
  # ---- Run ODE ----
  yini <- c(c(X0), c(bet0), xP0, yA0)
  times <- seq(0, 300, 0.1)
  parameters <- list(
    rP = rP, rA = rA, alpha = alpha, w = w, G = G, c = c, kP = kP, kA = kA,
    xm = xm, ym = ym, sd_m = sd_m, sd_kP = sd_kP, sd_kA = sd_kA,
    sd_cP = sd_cP, sd_cA = sd_cA, m = m
  )

  R.utils::withTimeout(
    solution <- deSolve::ode(y = yini, times = times, func = lotka, parms = parameters),
    timeout = timeout, onTimeout = "error"
  )

  # ---- Process solution ----
  solution <- as.data.frame(solution)
  X <- solution %>%
    select((2:(M + N + 1))) %>%
    set_names(c(paste0("P", 1:M), paste0("A", 1:N)))

  ForEffMatA <- solution %>%
    select((M + N + 2):(M * N + M + N + 1)) %>%
    set_names(c(paste0("A", rep(1:N, each = M))))

  x_trait <- solution %>%
    select((M * N + M + N + 2):(M * N + 2 * M + N + 1)) %>%
    set_names(c(paste0("P", 1:M)))

  y_trait <- solution %>%
    select((M * N + 2 * M + N + 2):(M * N + 2 * M + 2 * N + 1)) %>%
    set_names(c(paste0("A", 1:N)))

  # ---- Identify survivors ----
  survival <- X %>%
    tail(5) %>%
    colSums() %>%
    keep(~ .x > 1e-5) %>%
    names()

  P_sur <- survival[str_detect(survival, "^P")]
  A_sur <- survival[str_detect(survival, "^A")]

  # ---- Extract surviving densities and traits ----
  P_density <- as.matrix(X[, colnames(X) %in% P_sur, drop = FALSE])
  A_density <- as.matrix(X[, colnames(X) %in% A_sur, drop = FALSE])
  x_trait_sur <- as.matrix(x_trait[, colnames(x_trait) %in% P_sur, drop = FALSE])
  y_trait_sur <- as.matrix(y_trait[, colnames(y_trait) %in% A_sur, drop = FALSE])

  # ---- Network structure ----
  X<-as.matrix(solution[,2:(M+N+1)])
  ForEffMatA<-as.matrix(solution[,(M+N+2):(M*N+M+N+1)])
  net<-matrix(tail(ForEffMatA,1),M,N)*(tail(X,1)[1:M]%*%t(tail(X,1)[(M+1):(M+N)]))
  net[net < 1e-4] <- 0
  if(suffix==""){
    struct <- bipartite::networklevel(net, index = c("modularity","weighted nestedness","weighted connectance", "Shannon diversity"))
    names(struct) <- paste0(names(struct), suffix)

  }else{
    struct <- bipartite::networklevel(net, index = c("modularity","weighted nestedness","weighted connectance", "Shannon diversity"))
    names(struct) <- paste0(names(struct), suffix)
  }


  # Final plant competition coefficient
  disMatPf<-matrix(unlist(tail(x_trait,1)),M,M)-matrix(unlist(tail(x_trait,1)),M,M,byrow = T)
  cPf<-1*exp(-(disMatPf^2)/(2*(sd_cP)^2))
  mean_cPf <- mean(cPf[-seq(1, M*M, by = M + 1)])
  sd_cPf <- sd(cPf[-seq(1, M*M, by = M + 1)])

  # final animal competition coefficient
  disMatAf<-matrix(unlist(tail(y_trait,1)),N,N, byrow = T)-matrix(unlist(tail(y_trait,1)),N,N)
  cAf<-1*exp(-(disMatAf^2)/(2*(sd_cA)^2))
  mean_cAf <- mean(cAf[-seq(1, N*N, by = N + 1)])
  sd_cAf<-sd(cAf[-seq(1, N*N, by = N + 1)])


  # final plant and animal benefits
  disMatf=c(matrix(unlist(tail(x_trait,1)),M,1))-matrix(unlist(tail(y_trait,1)),M,N, byrow = T)
  sigmaf<-c*exp(-(disMatf^2)/(2*(sd_m)^2))
  mean_sigmaf<-mean(sigmaf)
  sd_sigmaf<-sd(sigmaf)

  #Final plant carrying capacity
  KPf<-kP*exp(-((unlist(tail(x_trait,1))-xm)^2)/(2*(sd_kP)^2))
  mean_KPf<-mean(KPf)
  sd_KPf<-sd(KPf)

  # Final animal carrying capacity
  KAf<-kA*exp(-((unlist(tail(y_trait,1))-ym)^2)/(2*(sd_kA)^2))
  mean_KAf<-mean(KAf)
  sd_KAf<-sd(KAf)

  # Steady state
  y_final<-unlist(c(tail(X,1),tail(ForEffMatA,1),tail(x_trait,1),tail(y_trait,1)))

  resilience <- compute_resilience(y_final, parameters)

  # ---- Summary stats ----
  stats <- c(
    persis_P = length(P_sur) / M,
    persis_A = length(A_sur) / N,
    persis_T = (length(P_sur)+length(A_sur)) / (M+N),
    max_x = max(tail(x_trait_sur, 1)), max_y = max(tail(y_trait_sur, 1)),
    mean_x = mean(tail(x_trait_sur, 1)), mean_y = mean(tail(y_trait_sur, 1)),
    sd_x = sd(c(tail(x_trait_sur, 1))), sd_y = sd(c(tail(y_trait_sur, 1))),
    rao_P = rao_quadratic_entropy(tail(P_density,1), tail(x_trait_sur,1)),
    rao_A = rao_quadratic_entropy(tail(A_density,1), tail(y_trait_sur,1)),
    even_P = pielou_evenness(tail(P_density,1)),
    even_A = pielou_evenness(tail(A_density,1)),
    mean_P = mean(c(tail(P_density, 1))), mean_A = mean(c(tail(A_density, 1))),
    sd_P = sd(c(tail(P_density, 1))), sd_A = sd(c(tail(A_density, 1))),
    sum_P = sum(c(tail(P_density, 1))), sum_A = sum(c(tail(A_density, 1))),
    w_sum_P = sum(c(tail(P_density, 1)) * c(tail(x_trait_sur, 1))),
    w_sum_A = sum(c(tail(A_density, 1)) * c(tail(y_trait_sur, 1))),
    min_P = min(c(tail(P_density, 1))), min_A = min(c(tail(A_density, 1))),
    mean_cPf = mean_cPf,
    mean_cAf = mean_cAf, mean_sigmaf = mean_sigmaf,
    sd_cPf = sd_cPf ,sd_cAf = sd_cAf, sd_sigmaf = sd_sigmaf,
    mean_KPf = mean_KPf, mean_KAf = mean_KAf, sd_KPf = sd_KPf, sd_KAf = sd_KAf,
    cv_P = mean(apply(tail(P_density,100),2,sd)/colMeans(tail(P_density,100)),na.rm=T),
    cv_A = mean(apply(tail(A_density,100),2,sd)/colMeans(tail(A_density,100)),na.rm=T),
    cv_x = mean(apply(tail(x_trait_sur,100),2,sd)/colMeans(tail(x_trait_sur,100)),na.rm=T),
    cv_y = mean(apply(tail(y_trait_sur,100),2,sd)/colMeans(tail(y_trait_sur,100)),na.rm=T),
    resilience = resilience
  )
  names(stats) <- paste0(names(stats), suffix)

  return(c(stats, struct))
}


model<-function(M,N){

  tryCatch(
    expr = {
      # Initial trait
      xP0=runif(M,1,2) # plant competitive trait
      yA0=runif(N,1,2) # animal competitive trait

      # optimal trait for maximum carrying capacity
      xm=3 # plant
      ym=3 # animal

      #mutational rate
      m=0.0005

      # standard deviation for mutualistic interaction
      km=sample(seq(-2,2,0.25),1)
      sd_m=exp(km)
      # standard deviation for carrying capacity
      kk=sample(seq(-2,2,0.25),1)
      sd_kP=exp(kk)
      sd_kA=exp(kk)

      # standard deviation for within guild competition kernel
      kc=sample(seq(-2,2,0.25),1)
      sd_cP=exp(kc)
      sd_cA=exp(kc)

      # maximum carrying capacity
      kP=10
      kA=10

      # Intrinsic growth rate
      rP=matrix(runif(M,0.1,0.2), nr=M)
      rA=matrix(runif(N,0.1,0.2), nr=N)

      # floral resource production rate
      alpha=10*matrix(runif(M,1,1.2), nr=M)

      # floral resource decay rate
      w=1*matrix(runif(M,1,1.2), nr=M)

      # Benefit scaling parameter

      c=1

      #Initial foraging effort
      bet0<-matrix(1/M,M,N)

      #Initial densities
      XP0<-matrix(runif(M,1,1.2), nr=M)
      XA0<-matrix(runif(N,1,1.2), nr=N)

      #Animal adaptation rate
      G=2*runif(N,1,1.2)
      # combined initial density
      X0=rbind(XP0,XA0)

      # Simulation of the model equation
      # Evo-Beh
      value1 <- run_scenario(M, N, lotka, X0, bet0, xP0, yA0,
                             rP, rA, alpha, w, c,
                             kP, kA, xm, ym,
                             sd_m, sd_kP, sd_kA, sd_cP, sd_cA,
                             G = G,  m = m, suffix = "")

      # Evo
      value2 <- run_scenario(M, N, lotka, X0, bet0, xP0, yA0,
                             rP, rA, alpha, w, c,
                             kP, kA, xm, ym,
                             sd_m, sd_kP, sd_kA, sd_cP, sd_cA,
                             G = G * 0, m = m, suffix = "1")

      # Beh
      value3 <- run_scenario(M, N, lotka, X0, bet0, xP0, yA0,
                             rP, rA, alpha, w, c,
                             kP, kA, xm, ym,
                             sd_m, sd_kP, sd_kA, sd_cP, sd_cA,
                             G = G,  m = 0, suffix = "2")

      # Baseline
      value4 <- run_scenario(M, N, lotka, X0, bet0, xP0, yA0,
                             rP, rA, alpha, w, c,
                             kP, kA, xm, ym,
                             sd_m, sd_kP, sd_kA, sd_cP, sd_cA,
                             G = G * 0, m = 0, suffix = "3")



      # Initial plant competition coefficient
      disMatP0<-matrix(xP0,M,M)-matrix(xP0,M,M,byrow = T)
      cP0<-1*exp(-(disMatP0^2)/(2*(sd_cP)^2))
      mean_cP0 <- mean(cP0[-seq(1, M*M, by = M + 1)])
      sd_cP0<-sd(cP0[-seq(1, M*M, by = M + 1)])


      # Initial animal competition coefficient
      disMatA0<-matrix(yA0,N,N, byrow = T)-matrix(yA0,N,N)
      cA0<-1*exp(-(disMatA0^2)/(2*(sd_cP)^2))
      mean_cA0 <- mean(cA0[-seq(1, N*N, by = N + 1)])
      sd_cA0<-sd(cA0[-seq(1, N*N, by = N + 1)])


      #Initial plant and animal benefits
      disMat0=c(matrix(xP0,M,1))-matrix(yA0,M,N, byrow = T)
      sigma0<-c*exp(-(disMat0^2)/(2*(sd_m)^2))
      mean_sigma0<-mean(sigma0)
      sd_sigma0<-sd(sigma0)

      # Initial plant carrying capacity
      KP0<-kP*exp(-((xP0-xm)^2)/(2*(sd_kP)^2))
      mean_KP0<-mean(KP0)
      sd_KP0<-sd(KP0)
      # Initial animal carrying capacity
      KA0<-kA*exp(-((yA0-ym)^2)/(2*(sd_kA)^2))
      mean_KA0<-mean(KA0)
      sd_KA0<-sd(KA0)


      # Combine all
      k <- c(km = km, kk = kk, kc = kc, P = M, A= N, alpha = mean(alpha),
             mean_x0=mean(xP0), mean_y0=mean(yA0),
             sd_x0 = sd(xP0), sd_y0 = sd(yA0), mean_cP0 = mean_cP0,
             mean_cA0 = mean_cA0, mean_sigma0 = mean_sigma0,
             sd_cP0 = sd_cP0 ,sd_cA0 = sd_cA0, sd_sigma0 = sd_sigma0,
             mean_KP0 = mean_KP0, sd_KP0 = sd_KP0, mean_KA0 = mean_KA0,
             sd_KA0 = sd_KA0)
      values <- c(k, value1, value2, value3, value4)
      return(values)

    },
    error = function(e){
      return(NA)
    }
  )

}

# Attached packages
# - dplyr       (v.1.4)
# - purrr       (v1.2.2)
# - stringr     (v1.6.0)
# - tidyr       (v1.3.1)
# - numDeriv    (v2016.8.1.1)

# Define the parallel loop
clusterEvalQ(cl,
             {library(dplyr)
               library(purrr)
               library(stringr)
               library(tidyr)
               library(numDeriv)})
clusterExport(cl,c("lotka",
                   "run_scenario",
                   "rao_quadratic_entropy",
                   "pielou_evenness",
                   "lotka_deriv",
                   "compute_resilience"))
system.time(trait_data <- foreach(M=rep(5:10,1),
                                  .inorder = FALSE) %:%
              foreach(N=rep(5:10,1)) %dopar% {
                result<-model(M,N)
              })

 # Stop cluster
stopCluster(cl)

trait_matrix<-data.frame()
for(i in 1:length(trait_data)){
  for(j in 1:length(trait_data[[i]])){
    if(length(trait_data[[i]][[j]])>1){
      trait_matrix<-rbind(trait_matrix,trait_data[[i]][[j]])
      names(trait_matrix)<-names(trait_data[[i]][[j]])
    }
  }
}

