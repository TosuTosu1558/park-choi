#Baker(1992)に基づいてexactにmlを求める
exact_ml_cal <-function(data){

	#パラメータ推定
	#Sayan Ghosh (2020) remark5.1
	A<-matrix(y_obs,I,J)
	A_D<-diag(diag(A),I,J)

	alpha<-solve(A_D)%*%y_mis_row
	ipshi_al<-10 #誤差

	while(abs(ipshi_al)>0.0000000001){
		alpha_old <- alpha
		alpha_new <- alpha_old + solve(A_D)%*%(y_mis_row-A%*%alpha_old)
		ipshi_al<-sum(alpha_new-alpha_old)
		alpha<-alpha_new
	}
	
	#Sayan Ghosh (2020)と一致
	#print(alpha)	

	beta<-solve(A_D)%*%y_mis_col
	ipshi_be<-10

	while(abs(ipshi_be)>0.0000000001){
		beta_old <- beta
		beta_new <- beta_old + solve(A_D)%*%(y_mis_col-t(A)%*%beta_old)
		ipshi_be<-sum(beta_new-beta_old)
		beta<-beta_new
	}

	#解が負になったら0を代入
	for (i in 1:I){
		if(alpha[i]<0){
			alpha[i]<-0
		}
		if(beta[i]<0){
			beta[i]<-0
		}
	}

	#パラメータ修正
	alpha[1]<-y_s21/sum(y_obs[1:I])
	g<-(sum(y_obs[1:I])*y_s22)/(y_s21*y_mis_col[1])

	#ML計算
	#欠測カテゴリー 11,12,21,22
	ml_11<-sum(y_obs[1:I])*(y_obs[1:I]+y_mis_row)/(sum(y_obs[1:I])+y_s21)
	ml_11<-append(ml_11,y_obs[3:4])

	A<-matrix(as.matrix(ml_11),2,2)
	A_D<-diag(diag(A),I,J)
	beta<-solve(A_D)%*%y_mis_col
	ipshi_be<-10

	while(abs(ipshi_be)>0.0000000001){
		beta_old <- beta
		beta_new <- beta_old + solve(A_D)%*%(y_mis_col-t(A)%*%beta_old)
		ipshi_be<-sum(beta_new-beta_old)
		beta<-beta_new
	}
	ml_12<- ml_11*(c(rep(beta,2)))
	

	ml_21<-ml_11*(c(rep(alpha[1],2),rep(alpha[2],2)))
	ml_22<-ml_11*(c(rep(beta,2)))*(c(rep(alpha[1],2),rep(alpha[2],2)))*g
	exact_ml<<-c(ml_11,ml_12,ml_21,ml_22)

	#尤度比統計検定量
	#Baker(1992)の結果と一致
	G<-2*(sum(y_obs*log(ml_11/y_obs))+
		sum(y_mis_col*log(c(ml_12[1]+ml_12[2],ml_12[3]+ml_12[4])/y_mis_col))+
		sum(y_mis_row*log(c(ml_21[1]+ml_21[3],ml_21[2]+ml_21[4])/y_mis_row))+
		y_s22*log(sum(ml_22)/y_s22))
	#print(G)
	
	#表出力
	ml_matrix_11<-matrix(as.matrix(exact_ml[1:4]),2,2,byrow=T)
	ml_matrix_12<-matrix(as.matrix(exact_ml[5:8]),2,2,byrow=T)
	ml_matrix_up<-cbind(ml_matrix_11,ml_matrix_12)
	ml_matrix_21<-matrix(as.matrix(exact_ml[9:12]),2,2,byrow=T)
	ml_matrix_22<-matrix(as.matrix(exact_ml[13:16]),2,2,byrow=T)
	ml_matrix_low<-cbind(ml_matrix_21,ml_matrix_22)
	ml_matrix<-rbind(ml_matrix_up,ml_matrix_low)
	#print(ml_matrix)

	rownames(ml_matrix)<-c("X_1=1","X_1=2","X_1=1","X_1=2")
	colnames(ml_matrix)<-c("X_2=1","X_2=2","X_2=1","X_2=2")
}





