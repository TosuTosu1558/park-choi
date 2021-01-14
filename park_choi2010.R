data <- read.table("Relationship between smoking and newborn weight.csv",
header=FALSE,sep=",")

data<-matrix(as.matrix(data),nrow(data),ncol(data))
y_obs<- matrix(t(data[-nrow(data),-ncol(data)]),
	1,length(data[-nrow(data),-ncol(data)]))
y_obs_matrix<-matrix(y_obs,2,2,byrow=T)
y_mis_col<-c(data[-nrow(data),ncol(data)])
y_mis_row<-c(data[nrow(data),-ncol(data)])
y_mis<-c(y_mis_col,y_mis_row,data[nrow(data),ncol(data)])
y_s12<-sum(as.matrix(data[1:nrow(data),ncol(data)]))-data[nrow(data),ncol(data)]
y_s21<-sum(as.matrix(data[nrow(data),1:ncol(data)]))-data[nrow(data),ncol(data)]
y_s22<-y_mis[length(y_mis)]
I<-nrow(data)-1
J<-ncol(data)-1
y_sum<-c(rep(sum(y_obs),I*2),rep(y_s12,I*2),rep(y_s21,I*2),rep(y_s22,I*2))
p<-(nrow(data)+1)*(ncol(data)+1)
N<-sum(data)

G_matrix<-c()

#Œv‰æs—ñ
Z<-matrix(c(1,1,1,1,1,1,1,1,1,
		1,1,0,1,1,1,0,0,1,
		1,0,1,1,1,0,1,0,1,
		1,0,0,1,1,0,0,0,1,
		1,1,1,1,0,1,0,1,0,
		1,1,0,1,0,1,0,0,0,
		1,0,1,1,0,0,0,0,0,
		1,0,0,1,0,0,0,0,0,
		1,1,1,0,1,0,1,1,0,
		1,1,0,0,1,0,0,0,0,
		1,0,1,0,1,0,1,0,0,
		1,0,0,0,1,0,0,0,0,
		1,1,1,0,0,0,0,1,0,
		1,1,0,0,0,0,0,0,0,
		1,0,1,0,0,0,0,0,0,
		1,0,0,0,0,0,0,0,0),
		16,9,byrow=T)


#•ªÍ
source("exact_ml_cal.R")
exact_ml_cal(data)
source("prior_cal.R")
prior_cal(data,exact_ml)
source("newton_raphson.R")
source("solve_matrix.R")
source("em.R")
em(data)


#GŒvZ
for(i in 1:6){
	ml_s12<-c(rep(c(sum(ml_matrix[i,(K+1):(3/2*K)]),
		sum(ml_matrix[i,(3/2*K+1):(2*K)])),each=2))
	ml_s21<-c(rep(c(sum(ml_matrix[i,seq(2*K+1,3*K,I)]),
		sum(ml_matrix[i,seq(2*K+2,3*K,I)])),2))
	ml_s22<-sum(ml_matrix[i,(3*K+1):(4*K)])

	G<- sum(y_obs*(log(ml_matrix[i,1:K]/y_obs))) + 
		sum(c(rep(y_mis_col,each=2))*log(ml_s12/c(rep(y_mis_col,each=2))))+
		sum(c(rep(y_mis_row,2))*log(ml_s21/c(rep(y_mis_row,2)))) +
		y_s22*log(ml_s22/y_s22)
	G_matrix<-append(G_matrix,G)
}

print(G_matrix)
print(ml_matrix)
print(Beta_matrix)










