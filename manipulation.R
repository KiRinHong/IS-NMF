
sinebell <- function(wl, ovlp){
  tmp = rep(0, wl)
  tmp[1:ovlp] = sin(pi*(1:ovlp - 0.5)/(2*ovlp))
  tmp[(ovlp+1):(wl-ovlp)] = 1
  tmp[(wl-ovlp+1):wl] = sin(pi*(wl-(wl-ovlp+1):wl + 0.5)/(2*ovlp))
  tmp
}
make_frames <- function(x, wl, ovlp){
  W = length(wl)
  NS = length(x)
  n_frames = ceiling((NS + ovlp)/(W - ovlp))
  Tpad = ovlp + n_frames*(W - ovlp)
  xpad = c(rep(0, ovlp), x, rep(0, Tpad-NS-ovlp))
  frames_index = 1 + 0:(n_frames-1)*(W-ovlp)
  Fx = matrix(0, W, n_frames)
  
  for (i in 1:n_frames) {
    Fx[,i] = t(xpad[frames_index[i]+0:(W-1)]*wl)
  }
  list(Fx = Fx, xpad = xpad)
}
stft <- function(x, wl, ovlp){
  if(length(wl) == 1){
    wl = sinebell(wl, ovlp)
  }
  W = length(wl)
  NS = length(x) # number of samples
  Fx = make_frames(x, wl, ovlp)$Fx
  STFT = apply(Fx, 2, fft) # treats the columns as vectors and returns the Fourier transform of each column
  if(W %% 2 == 0){
    STFT = STFT[1:(W/2+1),]
  }else{
    STFT = STFT[1:(W/2+1/2),]
  }
  STFT
}
overlap_add <- function(Fx, wl, ovlp){
  W = length(wl)
  n_frames = dim(Fx)[2]
  Tpad = ovlp + n_frames*(W - ovlp)
  xpad = rep(0, Tpad)
  frames_index = 1 + 0:(n_frames-1)*(W-ovlp)
  xpad[frames_index[1]:(frames_index[1]+W-1)] = Fx[,1] * wl
  for (i in 2:n_frames) {
    xpad[frames_index[i]:(frames_index[i]+W-1)] = xpad[frames_index[i]:(frames_index[i]+W-1)] + Fx[,i] * wl
  }
  xpad
}

man_istft <- function(STFT, wl, ovlp){
  if(length(wl) == 1){
    wl = sinebell(wl, ovlp)
  }
  W = length(wl)
  if(W %% 2 == 0){
    STFT = rbind2(STFT,Conj(STFT[(W/2):2,]))
  }else{
    STFT = rbind2(STFT,Conj(STFT[((W+1)/2):2,]))
  }
  ISTFT = Re(apply(STFT, 2, function(x) fft(x, inverse = TRUE)/length(x)))
  xpad = overlap_add(ISTFT, wl, ovlp)
  xpad
}
