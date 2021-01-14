my_solve <- function(A) {

#入力データチェック
len<-ncol(A)	
det_A <- det(A) # 行列式
if(det_A == 0) return ("正則行列じゃないよ")

	
#2×2行列の場合、余因子行列の行列式が求められないため
if(len == 2) return (matrix(1/det_A * c(A[2, 2], -A[1, 2], -A[2, 1], A[1, 1]), 2, 2, T))

#出力用逆行列生成
A_ <- matrix(0, len, len)
    
#余因子行列作成
for(i in 1 : len) {
	for(j in 1 : len) {
      	A_[j, i] <- ((-1)^(i + j)) * det(A[-i, -j])
      }
  }

return (A_ / det_A)
}