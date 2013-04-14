var('d,y,x')
def kronecker(l,m,n):
    s=SymmetricFunctionAlgebra(QQ,basis="schur")
    sl=s(l)
    sm=s(m)
    sn=s(n)
    k=sn.scalar(sl.itensor(sm))
    return k

def region(l1,l2,mu1,mu2,mu3,nu1,nu2):

    print "lambda: [",l1,l2,"]"
    print "mu:     [",mu1,mu2,mu3,"]"
    print "nu:     [",nu1,nu2,"]"
    
    
    if (mu1 < mu2) or (mu2 < mu3):
        print("mu no es una partición")
    elif not(mu1+mu2+mu3 == nu1+nu2):
        print("sumas inconsistentes: mu=",mu1+mu2+mu3, " nu=", nu1+nu2)
    elif not(mu1+mu2+mu3 == l1+l2):
        print "sumas inconsistentes"        
    else:
        print "coeff kronecker= ", kronecker( [l1,l2],[mu1,mu2,mu3],[nu1,nu2])
        
        nu1b = nu1 +1
        
        s= min(mu3,nu1b/3)
        q = min(nu1b+l1-mu1,nu1b+l2-mu2)
        
        xmax = 10
        xmin=0
        ymax=mu2+1
        ymin=mu3-1
        ymin=0
        
        #ineqs = [d>=0, y>=0, y<= t1 ,d+y-mu2 <= 0,d+y-mu3 >= 0,2*y+d<=l2,  nu1+l2-mu1<= 3*y+d, 3*y+d <= t3, 4*y+d <= nu1+mu3]
        
        # Desigualdades de nu barra bajo la transformaci'on
    
        
        
        q2=min(nu1,nu1+l1-mu1,nu1+l2-mu2)        
        
        ineqs = [ x >=0 , y >=x,  x <= s, y>=mu3, y<= mu2, x+y<=l2, 2*x+y>=nu1b+l2-mu1]+[3*x+y <= nu1b+mu3, 2*x+y <= q]
        G=region_plot(ineqs, (x,xmin,xmax), (y,ymin,ymax),incol="lightgray")
        
        ineqs = [x>=0, x<=mu3,y<=mu2,y>=mu3,x+y<=l2,x+2*y>=nu1+l2-mu1 ]+[x+2*y<=q2]+[2*x+y<nu1b+l2-mu1]
        G+=region_plot(ineqs, (x,xmin,xmax), (y,ymin,ymax),incol="orange")
        

        eqs =  [ x ==0 ,  y == mu3, y== mu2, x+y==l2]
        for eq in eqs:
            G+=implicit_plot(eq, (x,xmin,xmax), (y,ymin,ymax),color="black")
            
        eqs =  [  y ==x,  x == s,  2*x+y==nu1b+l2-mu1]
        for eq in eqs:
            G+=implicit_plot(eq, (x,xmin,xmax), (y,ymin,ymax),color="cyan")
        eqs =  [ 2*x+y==nu1b+l2-mu1,2*x+y == q]
        for eq in eqs:
            G+=implicit_plot(eq, (x,xmin,xmax), (y,ymin,ymax),color="blue")
        
        eqs =  [ 2*x+y==nu1+l2-mu1]
        for eq in eqs:
            G+=implicit_plot(eq, (x,xmin,xmax), (y,ymin,ymax),color="blue",linestyle="--")
                        
        eqs =  [ 3*x+y == nu1b+mu3]
        for eq in eqs:
            G+=implicit_plot(eq, (x,xmin,xmax), (y,ymin,ymax),color="cyan")            
            
        eqs = [  x==mu3,x+2*y==nu1+l2-mu1, x+2*y==q2]
        for eq in eqs:
            G+=implicit_plot(eq, (x,xmin,xmax), (y,ymin,ymax),color="red")            
            
        puntos=[(i,j) for i in range(xmin,xmax+1) for j in range(ymin,ymax+1)]
        #for p in puntos:
        #    G+=point(p,size=40)
        G+=list_plot(puntos)        
        rangox = range(xmin,xmax+1)
        rangoy = range(ymin,ymax+1)
        for lx in rangox:
            G+=line([(lx,0),(lx,ymax)],alpha=0.25,color="black")
        for ly in rangoy:
            G+=line([(0,ly),(xmax,ly)],alpha=0.25,color="black")
        return G

G=region(
l1=15,
l2=5,
mu1=12,
mu2=6,
mu3=2,
nu1=11,
nu2=9)
show(G)
