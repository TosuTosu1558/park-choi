source("newton_raphson.R")
source("solve_matrix.R")

em<-function(data){
	ml_matrix<-c()
	Beta_matrix<-c()
	ml_matrix<-c()
	iteration<-0

	for (i in 1:6){

		#初期値
		#MAR下のMLE
		m_t_11<-y_obs
		m_t_12<-(sum(y_mis_col)/sum(y_obs))*y_obs
		m_t_21<-(sum(y_mis_row)/sum(y_obs))*y_obs
		m_t_22<-(sum(y_s22)/sum(y_obs))*y_obs
		m_t<-c(m_t_11,m_t_12,m_t_21,m_t_22)
				
		likelihood<-sum(dmultinom(m_t,prob=m_t/N,log=T))
		dif<-10


		while(dif>0.00000003){

			#K=4

			#m_tの各欠測カテゴリの和
			m_t_s12<-c(rep(c(sum(m_t[(K+1):(3/2*K)]),sum(m_t[(3/2*K+1):(2*K)])),each=2))
			m_t_s21<-c(rep(c(sum(m_t[seq(2*K+1,3*K,I)]),sum(m_t[seq(2*K+2,3*K,I)])),2))
			m_t_s22<-sum(m_t[(3*K+1):(4*K)])


			#deltaの各欠測カテゴリの和
			delta_s12<-c(rep(c(sum(delta_matrix[i,(K+1):(3/2*K)]),
				sum(delta_matrix[(3/2*K+1):(2*K)])),each=2))
			delta_s21<-c(rep(c(sum(delta_matrix[i,seq(2*K+1,3*K,I)]),
				sum(delta_matrix[i,seq(2*K+2,3*K,I)])),2))
			delta_s22<-c(rep(sum(delta_matrix[i,(3*K+1):(4*K)]),4))

			#疑似観測値計算
			y_pes_11<-(sum(y_obs)/(sum(y_obs)+sum(delta_matrix[i,1:K])))*
				(y_obs+delta_matrix[i,1:K])
			y_pes_12<-(c(rep(y_mis_col,each=2))*(m_t[(K+1):(2*K)]/m_t_s12)+ 
				delta_matrix[i,(K+1):(2*K)])*
				(c(rep(y_mis_col,each=2)) /(c(rep(y_mis_col,each=2)+delta_s12)))
			y_pes_21<- (c(rep(y_mis_row,2))*(m_t[(2*K+1):(3*K)]/m_t_s21)
				+delta_matrix[i,(2*K+1):(3*K)])*
				(c(rep(y_mis_row,2))/(c(rep(y_mis_row,2))+delta_s21))
			y_pes_22<-y_s22*(m_t[(3*K+1):(4*K)]/m_t_s22)*
				(y_s22/(y_s22+delta_s22))
			y_pes<-c(y_pes_11,y_pes_12,y_pes_21,y_pes_22)


			#MPE計算
			Beta<-newton_raphson(m_t,y_pes)
			m_t<-as.vector(exp(Z%*%Beta))
			

			#誤差計算
			likelihood_new<-sum(dmultinom(m_t,prob=m_t/N,log=T))
			dif<-abs(likelihood-likelihood_new)
			likelihood<-likelihood_new
			iteration<-iteration+1
		}
		
		#推定結果を格納
		ml_matrix<-append(ml_matrix,exp(Z%*%Beta))
		Beta_matrix<-append(Beta_matrix ,Beta)

	}
	
	Beta_matrix<<-matrix(as.matrix(Beta_matrix),6,9,byrow=T)
	ml_matrix<<-matrix(as.matrix(ml_matrix),6,16,byrow=T)

	
	
}






