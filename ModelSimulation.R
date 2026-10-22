# R version
# R 4.5.2

# Attached packages
# - deSolve           (v1.40)
# - tidyverse         (v2.0.0)
# - bipartite         (v2.24)

# Packages required for ODE and network metrics
library(deSolve)    # Numerical solution of ordinary differential equations (ODEs)
library(tidyverse)  # Data manipulation, visualization, and data wrangling
library(bipartite)  # Analysis and visualization of bipartite ecological networks


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
    disMatP=outer(xP,xP,"-")
    disMatA=t(outer(yA,yA,"-"))
    # cross guild dissimilarity matrix
    disMat=outer(xP,yA,"-")
    # competition kernel

    cP=1*exp(-(disMatP^2)/(2*(sd_cP)^2))
    cA=1*exp(-(disMatA^2)/(2*(sd_cA)^2))


    # benefit function
    sigma_P=c*exp(-(disMat^2)/(2*(sd_m)^2))
    sigma_A=c*exp(-(disMat^2)/(2*(sd_m)^2))

    # Population dynamics/ Ecological process (Lotka-Volterra)
    dXP<-(rP-(rP*(cP%*%XP))/KP+(((sigma_P*bet)%*%XA)*(alpha))/(w+((sigma_A*bet)%*%XA)))
    dXA<-(rA-(rA*(cA%*%XA))/KA+(t(sigma_A*bet)%*%((alpha*XP)/(w+((sigma_A*bet)%*%XA)))))
    dX<-(Xx[1:(M+N)]*rbind(dXP,dXA))
    # Adaptive foraging (replicator equation)
    dbet<-((bet%*%diag(G,ncol = N))*sweep((sigma_A*c(alpha)*c(XP))/c(w+((sigma_A*bet)%*%XA)),2,
                                          colSums((bet*sigma_A*c(alpha)*c(XP))/c(w+((sigma_A*bet)%*%XA))),"-"))


    #  Selection gradient of plant trait
    gx=c(rP*XP/KP)*rowSums(cP*(-(xP-xm)/sd_kP^2+disMatP/sd_cP^2))+(((c(alpha)*bet^2*disMat*sigma_P*sigma_A)/c(sd_m^2*(((bet*sigma_A)%*%XA)+w)^2))%*%XA^2)-(((c(alpha)*bet*disMat*sigma_P)/c(sd_m^2*(((bet*sigma_A)%*%XA)+w)))%*%XA^2)


    # Selection gradient of animal trait
    gy=colSums((cA*(matrix((-(yA-ym)/sd_kA^2),N,N,byrow = TRUE)+disMatA/sd_cA^2))%*%diag(c((rA*XA/KA)),ncol = N))+colSums((c(XP*alpha)*bet*disMat*sigma_A)/c(sd_m^2*(((bet*sigma_A)%*%XA)+w)))

    # canonical equation of adaptive dynamics (trait dynamics)
    dxP=m*XP*gx
    dyA=m*XA*gy

    return(list(c(dX,dbet,dxP,dyA),cP,cA,sigma_P,sigma_A))
  })
}

# A function to extract the output of the simulation
solution_fun<-function(solution){
  solution<-as.data.frame(solution)
  # assume solution is a tibble or data.frame
  X <- solution %>%
    select((2:(M+N+1))) %>%
    set_names(c(paste0("P", 1:M), paste0("A", 1:N)))

  ForEffMatA <- solution %>%
    select((M+N+2):(M*N+M+N+1)) %>%
    set_names(c(paste0("A", rep(1:N, each = M))))

  x_trait <- solution %>%
    select((M*N+M+N+2):(M*N+2*M+N+1)) %>%
    set_names(c(paste0("P", 1:M)))

  y_trait <- solution %>%
    select((M*N+2*M+N+2):(M*N+2*M+2*N+1)) %>%
    set_names(c(paste0("A", 1:N)))

  # Species that survive (density > 1e-5 in last 5 timesteps)
  survival <- X %>%
    tail(5) %>%
    colSums() %>%
    keep(~ .x > 1e-5) %>%
    names()

  # Split into plants and animals
  P_sur <- survival[str_detect(survival, "^P")]
  A_sur <- survival[str_detect(survival, "^A")]

  X<-as.matrix(solution[,2:(M+N+1)])
  ForEffMatA<-as.matrix(solution[,(M+N+2):(M*N+M+N+1)])
  net<-matrix(tail(ForEffMatA,1),M,N)*(tail(X,1)[1:M]%*%t(tail(X,1)[(M+1):(M+N)]))
  net[net<0.0001]<-0

  # Set species names
  rownames(net) <- paste0("P", 1:M)
  colnames(net) <- paste0("A", 1:N)

  # Keep only surviving plants and animals
  net <- net[P_sur, A_sur, drop = FALSE]

  return(list("X"=X,
              "ForEffMatA"=ForEffMatA,
              "net"=net,
              "x_trait" = x_trait,
              "y_trait" = y_trait))

}


# No. of plant
M=10
# No. of animal
N=10

# Seed set
set.seed(12)


{  # Initial trait values
  xP0=runif(M,1,2) # plant competitive trait
  yA0=runif(N,1,2) # animal competitive trait
  # mutational rate
  m=0.0005
  # optimal trait for maximum carrying capacity
  xm=3 # plant
  ym=3 # animal


  # standard deviation for mutualistic interaction

  sd_m=exp(1)
  # standard deviation for carrying capacity

  sd_kP=exp(1)
  sd_kA=exp(1)

  # standard deviation for within guild competition kernel

  sd_cP=exp(-2)
  sd_cA=exp(-2)


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
  bet0 = matrix(1/M,M,N)

  #Initial densities
  XP0<-matrix(runif(M,1,1.2), nr=M)
  XA0<-matrix(runif(N,1,1.2), nr=N)

  #Animal adaptation rate
  G=2*runif(N,1,1.2)
  # combined initial density
  X0=rbind(XP0,XA0)
  yini=c(c(X0),c(bet0),xP0,yA0)}



#### No adaptation ECOLOGICAL (BASELINE) ####
# Simulation of the model equation
{times=seq(0,300,1)
  parameters=list(alpha=alpha,w=w, G=G*0,c=c, m = 0)
  solution1<-ode(y=yini, times=times, func=lotka, parms=parameters)}
solu1<-solution_fun(solution1)

#### EVOLUTIONARY ####
{times=seq(0,300,1)
parameters=list(alpha=alpha,w=w, G=G*0,c=c, m = m)
solution2<-ode(y=yini, times=times, func=lotka, parms=parameters)
}
solu2<-solution_fun(solution2)

#### BEHAVIOURAL ####
{times=seq(0,300,1)
parameters=list(alpha=alpha,w=w, G=G,c=c, m = 0)
solution3<-ode(y=yini, times=times, func=lotka, parms=parameters)}
solu3<-solution_fun(solution3)

#### EVOLUTIONARY AND BEHAVIOURAL ####
{times=seq(0,300,1)
  parameters=list(alpha=alpha,w=w, G=G,c=c, m = m)
  solution4<-ode(y=yini, times=times, func=lotka, parms=parameters)}
solu4<-solution_fun(solution4)


scens  <- list(solu1, solu2, solu3, solu4)
titles <- c("Baseline", "Evo", "Beh", "Evo-Beh")

layout(
  matrix(1:24, nrow = 6, byrow = TRUE),
  heights = c(1.3, 1, 1, 1, 1, 1)
)

par(
  cex.lab = 1.3,
  cex.axis = 1.5,
  cex.main = 1.8
)
## ========= Row 1 : Networks =========
for(i in 1:4){
  par(mar = c(1,1,3,1))
  plotweb(
    scens[[i]]$net,
    main = titles[i]
  )
}

## ========= Row 2 : Plant densities =========
for(i in 1:4){
  par(mar = c(2,4,1,1))
  matplot(
    scens[[i]]$X[,1:M],
    type = "l", lwd = 2, col = 1:M,
    xaxt = "n",
    main = "",lty = "solid" ,
    ylab = if(i == 1) "Plant density" else "",
    xlab = ""
  )
}

## ========= Row 3 : Pollinator densities =========
for(i in 1:4){
  par(mar = c(2,4,1,1))
  matplot(
    scens[[i]]$X[,(M+1):(M+N)],
    type = "l", lwd = 2, col = 1:N,
    xaxt = "n",
    main = "",lty = "solid" ,
    ylab = if(i == 1) "Pollinator density" else "",
    xlab = ""
  )
}

## ========= Row 4 : x trait =========
for(i in 1:4){
  par(mar = c(2,4,1,1))
  matplot(
    scens[[i]]$x_trait,
    type = "l", lwd = 2, col = 1:M,
    xaxt = "n",
    main = "",lty = "solid" ,
    ylab = if(i == 1) "Trait x" else "",
    xlab = ""
  )
}

## ========= Row 5 : y trait =========
for(i in 1:4){
  par(mar = c(2,4,1,1))
  matplot(
    scens[[i]]$y_trait,
    type = "l", lwd = 2, col = 1:N,
    xaxt = "n",
    main = "",lty = "solid" ,
    ylab = if(i == 1) "Trait y" else "",
    xlab = ""
  )
}

## ========= Row 6 : Foraging effort =========
for(i in 1:4){
  par(mar = c(4,4,1,1))
  matplot(
    scens[[i]]$ForEffMatA,
    type = "l", lwd = 2, col = 1:N,
    main = "",lty = "solid" ,
    ylab = if(i == 1) "Foraging effort" else "",
    xlab = "Time"
  )
}



