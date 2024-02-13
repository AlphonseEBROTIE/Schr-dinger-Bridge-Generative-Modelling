
# Algo Schrodinger Bridge TS , implementation Alphonse
# Novembre 2023


######### fonctions utiles pour l'implementation de l'algorithme

# fonction kernel

# K=function(x){
#   return((1-(abs(x))^2)*(abs(x)<=1))
# }

K_h_x=function(x,h){
  K=(1-(abs(x/h))^2)*(abs(x/h)<=1)
  return((1/h)*K)
}



######## Fonction Fi : utile pour l'appoximation proposee par Pham et al.

Fi=function(t,xi,y_k,xi_plu_un,i,t_k_i){
  t_i_plus_1=i+1
  t=t_k_i
  a=-((abs(xi_plu_un-y_k))^2)/(2*(t_i_plus_1-t))
  b=((abs(xi_plu_un-xi))^2)/(2*(t_i_plus_1-i))
  return(exp(a+b))
}


######## Fonction a : l'estimateur du path-dependent drift de la solution 
######## au pb de  controle stochastique


# a_hat=function(t,y_k,i,x,t_k_i){ #xi=x[i]
#   
#   t_i_plus_1 = i+1
#   
#   
#   u=sapply(1:i,function(j)K_h_x(x[j]-X_m[1:M,j],h)) # estimation par kernel des differents paths
#   
#   S=((X_m[,i+1])-y_k)*Fi(t,X_m[,i],y_k,X_m[,i+1],i,t_k_i)*sapply(1:M,function(j)prod(u[j,]))
#   
#   V=Fi(t,X_m[,i],y_k,X_m[,i+1],i,t_k_i)*sapply(1:M,function(j)prod(u[j,]))
#   
#   
#   res= sum(S)/(sum(V)*(t_i_plus_1-t))
#   
#   if(res=='NaN'){return(0)}else{return(res)}
#   
# }


# version plus rapide de a_hat

a_hat = function(t, y_k, i, x, t_k_i){
  
  t_i_plus_1 = i + 1
  u = sapply(1:i, function(j) K_h_x(x[j] - X_m[1:M, j], h))

  delta_X = X_m[, i + 1] - y_k
  fi_values = Fi(t, X_m[, i], y_k, X_m[, i + 1], i, t_k_i)
  prod_u = apply(u, 1, prod)

  S = delta_X * fi_values * prod_u
  V = fi_values * prod_u

  res = sum(S) / (sum(V) * (t_i_plus_1 - t))

  if(is.nan(res)) { return(0) } else { return(res) }
  
}





SBTS=function(X_m,n_simu,N_pi,h){
  M=nrow(X_m)
  N=ncol(X_m)-1
  Matrice <- matrix(0, nrow = n_simu, ncol = N)
  time_start=Sys.time()

  for (r in 1:n_simu){

    x <- rep(0,N + 1)
    y <- rep(0,N_pi + 1)

    for(i in 1:N){

      y[1]=x[i]

      for(k in 1:N_pi){

        t_k_i <- i + ((k - 1) / N_pi)

        #schema d'euler

        y[k+1]=y[k]+(1/N_pi)*a_hat(t_k_i,y[k],i,x,t_k_i) + (1/sqrt(N_pi))*rnorm(1,0,1)

      }
      x[i+1]=y[N_pi+1]
    }

    #matrice qui recupere les differents paths generes

    Matrice[r, ] <- x[-1]

  }


  generated <- as.data.frame(Matrice)

  time_end=Sys.time()

  print(time_end-time_start)

  return(generated)

}




########################################################################################
#########################  EXPERIMENTATIONS NUMERIQUES  ################################




#modele : toy autoregressive model

b=0.7;sig1=0.1;sig2=sig3=0.05
beta1=beta2=-1
M=1000

X_0=0
X_t1= b + rnorm(M,0,sig1)
X_t2= (beta1*X_t1) + rnorm(M,0,sig2)
X_t3 = (beta2*X_t2) +sqrt(abs(X_t1)) + rnorm(M,0,sig3)


X_m=cbind(X_0,as.data.frame(X_t1),as.data.frame(X_t2),as.data.frame(X_t3))



# visualisation des donnees

#plot(X_m$X_t1, X_m$X_t2, pch=19, col="blue")




#### compute


n_simu=500
N_pi=500

#h=0.05 dans le papier

# estimateur du bandwidth :

#h=bw.SJ(c(X_0,X_t1,X_t2,X_t3))


#methode 2 : pour data en multivarie 

library(KernSmooth)
h=dpik(c(X_0,X_t1,X_t2,X_t3));h

X_m=X_m


X_simu=SBTS(X_m=X_m,n_simu=n_simu,N_pi=N_pi,h=h)

sum(X_simu$V3<0)
X=X_simu

# X=X_simu[-which(X_simu$V3<0),]
# X=X[-which(X$V3<1),]
# 
# #X=X[which( (X$V3 !=max(X$V3)) & X$V2 !=max(X$V2)),]

# Creer des nuages de points pour les simulations

plot(X$V1,X$V3,pch=19, col="red")


points(X_t1, X_t3, pch = 19, col = "blue")
#plot(X_t1,X_t2)

plot(density(X_t3),col="blue")

lines(density(X_simu$V3),col="red")

polygon(density(X$V3), col =rgb(0, 0, 1))
polygon(density(X_t3), col =  rgb(1, 0, 0))


### difference de correlation

cor(X$V1,X$V2)-cor(X_t1,X_t2)
cor(X$V1,X$V3)-cor(X_t1,X_t3)
cor(X$V2,X$V3)-cor(X_t2,X_t3)


### comparaison quantiles (q_5 , q_95)

quantiles <- quantile(X_t1, probs = c(0.05, 0.95))
quantiles

quantiles1 <- quantile(X$V1, probs = c(0.05, 0.95))
quantiles1

quantiles2 <- quantile(X_t2, probs = c(0.05, 0.95))
quantiles2

quantiles3 <- quantile(X$V2, probs = c(0.05, 0.95))
quantiles3


quantiles4 <- quantile(X_t3, probs = c(0.05, 0.95))
quantiles4

quantiles5 <- quantile(X$V3, probs = c(0.05, 0.95))
quantiles5
