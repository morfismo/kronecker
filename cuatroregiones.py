var('d,y,x')
def kronecker(l,m,n):
    s=SymmetricFunctionAlgebra(QQ,basis="schur")
    sl=s(l)
    sm=s(m)
    sn=s(n)
    k=sn.scalar(sl.itensor(sm))
    return k

def regionP1(l1,l2,mu1,mu2,mu3,nu1,nu2):
  
    
    if (mu1 < mu2) or (mu2 < mu3):
        print("mu no es una partición")
    elif not(mu1+mu2+mu3 == nu1+nu2):
        print("sumas inconsistentes: mu=",mu1+mu2+mu3, " nu=", nu1+nu2)
    elif not(mu1+mu2+mu3 == l1+l2):
        print "sumas inconsistentes"        
    else:
        #print "coeff kronecker= ", kronecker( [l1,l2],[mu1,mu2,mu3],[nu1,nu2])
        
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
        

        #REGION P1
        ineqs = [x>=0, x<=mu3,y<=mu2,y>=mu3,x+y<=l2,x+2*y>=nu1+l2-mu1 ]+[x+2*y<=q2]
        
        #REGION P2
        # g0 = max(0, l2-mu1)3m
        # g1 =min(l1-mu1, l2-mu2)
        #ineqs = [x>=g0, x <= g1 , y >= nu1/3, x+2*y <= nu1, x+y <= mu2,  x+y >= mu3, y-x >= nu1-l2, 2y-x>=nu1-mu3]
        
        
        
        G=region_plot(ineqs, (x,xmin,xmax), (y,ymin,ymax),incol="blue")
        
        eqs =  [ x ==0 ,  x == mu3, y== mu2, y==mu3, x+y==l2,x+2*y==nu1+l2-mu1 ,x+2*y==q2]
        for eq in eqs:
            G+=implicit_plot(eq, (x,xmin,xmax), (y,ymin,ymax),color="black")
            
        puntos=[(i,j) for i in range(xmin,xmax+1) for j in range(ymin,ymax+1)]
        #for p in puntos:
        #    G+=point(p,size=40)
        G=G+list_plot(puntos)        
        rangox = range(xmin,xmax+1)
        rangoy = range(ymin,ymax+1)
        for lx in rangox:
            G+=line([(lx,0),(lx,ymax)],alpha=0.25,color="black")
        for ly in rangoy:
            G+=line([(0,ly),(xmax,ly)],alpha=0.25,color="black")
        return G
        
def regionP2(l1,l2,mu1,mu2,mu3,nu1,nu2):


    
    
    if (mu1 < mu2) or (mu2 < mu3):
        print("mu no es una partición")
    elif not(mu1+mu2+mu3 == nu1+nu2):
        print("sumas inconsistentes: mu=",mu1+mu2+mu3, " nu=", nu1+nu2)
    elif not(mu1+mu2+mu3 == l1+l2):
        print "sumas inconsistentes"        
    else:
        #print "coeff kronecker= ", kronecker( [l1,l2],[mu1,mu2,mu3],[nu1,nu2])
        
        nu1b = nu1 +1
        
        s= min(mu3,nu1b/3)
        q = min(nu1b+l1-mu1,nu1b+l2-mu2)
        
        xmax = 10
        xmin=0
        ymax=max(nu1-l2,(nu1-mu3)/2)+1
        ymin=mu3-1
        ymin=0
        
        #ineqs = [d>=0, y>=0, y<= t1 ,d+y-mu2 <= 0,d+y-mu3 >= 0,2*y+d<=l2,  nu1+l2-mu1<= 3*y+d, 3*y+d <= t3, 4*y+d <= nu1+mu3]
        
        # Desigualdades de nu barra bajo la transformaci'on
    
        
        q2=min(nu1,nu1+l1-mu1,nu1+l2-mu2)        
        

        
        #REGION P2
        g0 = max(0, l2-mu1)
        g1 =min(l1-mu1, l2-mu2)
        ineqs = [x>=g0, x <= g1 , y >= nu1/3, x+2*y <= nu1, x+y <= mu2,  x+y >= mu3, y-x >= nu1-l2, 2*y-x>=nu1-mu3]
        
        #REGION P1 cap P2
        #b1 = min(mu2,nu1/2)
        #b0 = max(nu-l2, (nu1-mu3)/2, nu1/3,mu3)
        
        G=region_plot(ineqs, (x,xmin,xmax), (y,ymin,ymax),incol="blue")
        
        eqs =  [ x==g0, x == g1 , y == nu1/3, x+2*y == nu1, x+y == mu2,  x+y == mu3]
        for eq in eqs:
            G+=implicit_plot(eq, (x,xmin,xmax), (y,ymin,ymax),color="black")

        eqs =  [ y-x== nu1-l2, 2*y-x==nu1-mu3]
        for eq in eqs:
            G+=implicit_plot(eq, (x,xmin,xmax), (y,ymin,ymax),color="blue")

            
        puntos=[(i,j) for i in range(xmin,xmax+1) for j in range(ymin,ymax+1)]
        #for p in puntos:
        #    G+=point(p,size=40)
        G=list_plot(puntos)+G
        rangox = range(xmin,xmax+1)
        rangoy = range(ymin,ymax+1)
        for lx in rangox:
            G+=line([(lx,0),(lx,ymax)],alpha=0.25,color="black")
        for ly in rangoy:
            G+=line([(0,ly),(xmax,ly)],alpha=0.25,color="black")
        return G

def regionPbarra(l1,l2,mu1,mu2,mu3,nu1,nu2):


    
    
    if (mu1 < mu2) or (mu2 < mu3):
        print("mu no es una partición")
    elif not(mu1+mu2+mu3 == nu1+nu2):
        print("sumas inconsistentes: mu=",mu1+mu2+mu3, " nu=", nu1+nu2)
    elif not(mu1+mu2+mu3 == l1+l2):
        print "sumas inconsistentes"        
    else:
        #print "coeff kronecker= ", kronecker( [l1,l2],[mu1,mu2,mu3],[nu1,nu2])
        
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
    
        q1=min(mu3,nu1b/3)
        q2=min(nu1,nu1+l1-mu1,nu1+l2-mu2)        
        q3 = min(nu1b+l1-mu1, nu1b+l2-mu2)
          
        #REGION Pbarra
        #b1 = min(mu2,nu1/2)
        #b0 = max(nu-l2, (nu1-mu3)/2, nu1/3,mu3)
        
        ineqs=[y<=q1, y>=0, x>=0, x+2*y<= l2, x+y<=mu2, x+y>=mu3, x+4*y<= (nu1b+mu3), x+3*y>=nu1b+l2-mu1, x+3*y<=q3]
        
        G=region_plot(ineqs, (x,xmin,xmax), (y,ymin,ymax),incol="blue")
        
        eqs =  [ y==q1]
        for eq in eqs:
            G+=implicit_plot(eq, (x,xmin,xmax), (y,ymin,ymax),color="gray")

        eqs =  [ x+2*y== l2]
        for eq in eqs:
            G+=implicit_plot(eq, (x,xmin,xmax), (y,ymin,ymax),color="black", linestyle="--")

        eqs =  [ x+4*y== (nu1b+mu3)]
        for eq in eqs:
            G+=implicit_plot(eq, (x,xmin,xmax), (y,ymin,ymax),color="black", linestyle="-.")

        eqs =  [  x+y==mu2, x+y==mu3 ]
        for eq in eqs:
            G+=implicit_plot(eq, (x,xmin,xmax), (y,ymin,ymax),color="red")
            
        eqs =  [ x+3*y==nu1b+l2-mu1, x+3*y==q3]
        for eq in eqs:
            G+=implicit_plot(eq, (x,xmin,xmax), (y,ymin,ymax),color="orange")
            
                                                
        puntos=[(i,j) for i in range(xmin,xmax+1) for j in range(ymin,ymax+1)]
        #for p in puntos:
        #    G+=point(p,size=40)
        G=list_plot(puntos)+G
        rangox = range(xmin,xmax+1)
        rangoy = range(ymin,ymax+1)
        for lx in rangox:
            G+=line([(lx,0),(lx,ymax)],alpha=0.25,color="black")
        for ly in rangoy:
            G+=line([(0,ly),(xmax,ly)],alpha=0.25,color="black")
        return G

        
def regiones(        l1,l2,mu1,mu2,mu3,nu1,nu2):
    
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
        pass
    
        #REGION P1 cap P2
        b1 = min(mu2,nu1/2)
        b0 = max(nu1-l2, (nu1-mu3)/2, nu1/3,mu3)        
        G1 = regionP1(l1,l2,mu1,mu2,mu3,nu1,nu2)
        G2 = regionP2(l1,l2,mu1,mu2,mu3,nu1,nu2)
        Gb = regionPbarra(l1,l2,mu1,mu2,mu3,nu1,nu2)

        show(G1)
        show(G2)
        
        print "interseccion:", floor(b1)-ceil(b0)+1
        
        show(Gb)