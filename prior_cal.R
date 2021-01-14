source("exact_ml_cal.R")
exact_ml_cal(data)

L<<-length(exact_ml)
K<<-L/4 #各欠測カテゴリーの長さ

#park and choi(2010)に基づいてHyper parameterを求める
#exact_ml_calを用いて求める

prior_cal<-function(data,exact_ml){

	ml_s11<-sum(exact_ml[1:K])
	ml_s12<-sum(exact_ml[K+1:K*2])
	ml_s21<-sum(exact_ml[(K*2+1):(K*3)])
	ml_s22<-sum(exact_ml[(K*3+1):(K*4)])
	ml_sum<-c(rep(ml_s11,K),rep(ml_s12,K),
		rep(ml_s21,K),rep(ml_s22,K))

	delta1 <<- (p/sum(data))/sum(y_obs) * (c(rep(y_obs,K))*y_sum)
	delta2 <<- (p/(sum(y_obs)*(sum(data)-sum(y_obs))))*y_sum*
		c(rep(0,4),rep(y_obs,3))
	delta3 <<- (p/(2*sum(exact_ml)))*ml_sum*
		c((2/ml_s11)*exact_ml[1:K],
		(exact_ml[(K+1):(2*K)]/ml_s12)+(1/(I*J)),
		(exact_ml[(2*K+1):(3*K)]/ml_s21)+(1/(I*J)),
		(exact_ml[(3*K+1):(4*K)]/ml_s22)+(1/(I*J)))
	delta4 <<-(p/(2*(sum(exact_ml)-ml_s11)))*ml_sum*
		c(rep(0,K),
		(exact_ml[(K+1):(2*K)]/ml_s12)+(1/(I*J)),
		(exact_ml[(2*K+1):(3*K)]/ml_s21)+(1/(I*J)),
		(exact_ml[(3*K+1):(4*K)]/ml_s22)+(1/(I*J)))
	delta5 <<-(p/3)*c(rep(0,K),rep(1/(I*J),3*K))

	delta_matrix<<-matrix(as.matrix(delta1,delta2,delta3,delta4,delta5),5,16)
	delta_matrix<<-rbind(c(rep(0,16)),delta_matrix)

}









