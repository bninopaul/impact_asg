#Niño Paul Batanay
#June 19, 2017

#Problem 9.6
def steepest_descent(Q, b, x0, eps=1e-10):
    norm_criterion = 1
    while norm_criterion > eps:
        Df = np.dot(x0,Q)-b
        Df_T = np.dot(Q,x0)-b
        alpha = np.dot(Df,Df_T)/np.dot(Df,np.dot(Q,Df_T))
        x1 = x0 - alpha*Df_T
        x0 = x1
        norm_criterion = np.sqrt(np.sum(Df**2))
    return x1

#Problem 9.7
def df_forward_diff(f, x, rerrf=1e-16):
    h = 2*np.sqrt(rerrf)
    Es = np.eye(len(x))
    Df = (f((x+h*Es).T) - f(x))/h
    return Df

