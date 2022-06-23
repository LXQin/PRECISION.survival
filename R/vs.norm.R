#' Variance stabilizing normalization

"vs.norm" <- function(train = NULL, test = NULL, ref.dis = NULL){
  if(!is.null(train) && !is.null(test)) stopifnot(nrow(train) == nrow(test))
  if(is.null(train)) stopifnot(!is.null(ref.dis))
  
  if(!is.null(train)){
    # vsn training
    train0 <- 2^train
    ref.dis <- vsn::vsn2(train0)
    train.vsn <- log2(exp(as.matrix(ref.dis)))
  } else train.vsn <- NULL
  
  if(is.null(test)) {
    test.fvsn <- NULL
  } else {
    test0 <- 2^test
    test.fvsn0 <- vsn::vsn2(test0, ref.dis)
    test.fvsn <- log2(exp(as.matrix(test.fvsn0)))
  }
  
  return(list("train.vsn" = train.vsn,
              "test.fvsn" = test.fvsn,
              "ref.dis" = ref.dis))
}
