var('d,y,x')
def kronecker(l,m,n):
	s=SymmetricFunctionAlgebra(QQ,basis="schur")
	sl=s(l)
	sm=s(m)
	sn=s(n)
	k=sn.scalar(sl.itensor(sm))
	return k

def regionP1(l1,l2,mu1,mu2,mu3,nu1,nu2):
	  
	nu1b = nu1 +1
		
	#s= min(mu3,nu1/3)
	#q = min(nu1+l1-mu1,nu1+l2-mu2)
		
	xmax = mu3+3
	xmin=0
	ymax=mu2+1
	ymin=mu3-1
	ymin=0
			
	q2=min(nu1,nu1+l1-mu1,nu1+l2-mu2)		 
		
	#REGION P1
	ineqs = [x>=0, x<=mu3,y<=mu2,y>=mu3,x+y<=l2,x+2*y>=nu1+l2-mu1 ]+[x+2*y<=q2]
		
	G=region_plot(ineqs, (x,xmin,xmax), (y,ymin,ymax),incol="blue")
		
	eqs =  [ x ==0 ,  x == mu3, y== mu2, y==mu3, x+y==l2,x+2*y==nu1+l2-mu1 ,x+2*y==q2]
	for eq in eqs:
		G+=implicit_plot(eq, (x,xmin,xmax), (y,ymin,ymax),color="black")
			
	rangox = range(xmin,xmax+1)
	rangoy = range(ymin,ymax+1)
	for lx in rangox:
		G+=line([(lx,0),(lx,ymax)],alpha=0.25,color="black")
	for ly in rangoy:
		G+=line([(0,ly),(xmax,ly)],alpha=0.25,color="black")
	return G
		
def regionP2(l1,l2,mu1,mu2,mu3,nu1,nu2):
	#REGION P2
	g0 = max(1, l2-mu1)
	g1 =min(l1-mu1, l2-mu2)
	
	xmax = mu2+1
	xmin=min(g0,g1)-1
	ymax=max(nu1-l2,(nu1-mu3)/2,mu2,nu1/2)+1
	ymin=mu3-1
	ymin=0
	
	q2=min(nu1,nu1+l1-mu1,nu1+l2-mu2)
	
	if g1 - g0 < 0 : print "La region P2 ES VACIA! (límite derecho es negativo)"

	

	ineqs = [x>=g0, x <= g1 , y >= nu1/3, x+2*y <= nu1, x+y <= mu2,  x+y >= mu3, y-x >= nu1-l2, 2*y-x>=nu1-mu3]

	G=region_plot(ineqs, (x,xmin,xmax), (y,ymin,ymax),incol="blue")
	
	eqs =  [ x==g0, x == g1 , x+2*y == nu1]
	for eq in eqs:
		G+=implicit_plot(eq, (x,xmin,xmax), (y,ymin,ymax),color="black")

	eqs =  [ x+y == mu2,  x+y == mu3]
	for eq in eqs:
		G+=implicit_plot(eq, (x,xmin,xmax), (y,ymin,ymax),color="darkgray")

	eqs =  [y == nu1/3 , y-x== nu1-l2, 2*y-x==nu1-mu3]
	for eq in eqs:
		G+=implicit_plot(eq, (x,xmin,xmax), (y,ymin,ymax),color="blue")

	rangox = range(xmin,xmax+1)
	rangoy = range(ymin,ymax+1)
	for lx in rangox:
		G+=line([(lx,ymin),(lx,ymax)],alpha=0.25,color="black")
	for ly in rangoy:
		G+=line([(xmin,ly),(xmax,ly)],alpha=0.25,color="black")
	return G

def regionPbarra(l1,l2,mu1,mu2,mu3,nu1,nu2):
    
	nu1b = nu1 +1
	
	s= min(mu3,nu1b/3)
	q = min(nu1b+l1-mu1,nu1b+l2-mu2)
	
	xmax = mu2+1
	xmin=0
	ymax=mu2+1
	ymin=0
	

	q1=min(mu3,nu1b/3)
	q2=min(nu1,nu1+l1-mu1,nu1+l2-mu2)		 
	q3 = min(nu1b+l1-mu1, nu1b+l2-mu2)
	  
	#REGION Pbarra
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
			
	rangox = range(xmin,xmax+1)
	rangoy = range(ymin,ymax+1)
	for lx in rangox:
		G+=line([(lx,ymin),(lx,ymax)],alpha=0.25,color="black")
	for ly in rangoy:
		G+=line([(xmin,ly),(xmax,ly)],alpha=0.25,color="black")
	return G

def regionPbarraT(l1,l2,mu1,mu2,mu3,nu1,nu2):
        # REGION TRANSFORMADA PARA P-BARRA
        
        
	nu1b = nu1 +1
	
	s= min(mu3,nu1b/3)
	q = min(nu1b+l1-mu1,nu1b+l2-mu2)
	
	xmax = mu3+3
	xmin=0
	ymax=mu2+1
	ymin=0
			

	

	q1=min(mu3,nu1b/3)		 
	q2 = min(nu1b+l1-mu1, nu1b+l2-mu2)
	  
	#REGION Pbarra-transformada
	ineqs=[y>=mu3, y<=mu2, x>=0, x<=q1, x+y<=l2, y>=x, 2*x+y <= q2,3*x+y<=nu1b+mu3,2*x+y>=nu1b+l2-mu1]


	G=region_plot(ineqs, (x,xmin,xmax), (y,ymin,ymax),incol="blue")
	
	eqs =  [y==mu3, y==mu2, x==0, x==q1, x+y==l2, 2*x+y == q2,3*x+y==nu1b+mu3]
	for eq in eqs:
		G+=implicit_plot(eq, (x,xmin,xmax), (y,ymin,ymax),color="darkgray")

	eqs =  [y==mu3, y==mu2, x==0, x==q1]
	for eq in eqs:
		G+=implicit_plot(eq, (x,xmin,xmax), (y,ymin,ymax),color="black")

	eqs =  [ y==x, 2*x+y==nu1b+l2-mu1]
	for eq in eqs:
		G+=implicit_plot(eq, (x,xmin,xmax), (y,ymin,ymax),color="blue")
			
	rangox = range(xmin,xmax+1)
	rangoy = range(ymin,ymax+1)
	for lx in rangox:
		G+=line([(lx,0),(lx,ymax)],alpha=0.25,color="black")
	for ly in rangoy:
		G+=line([(0,ly),(xmax,ly)],alpha=0.25,color="black")
	return G

		
def regiones(		 l1,l2,mu1,mu2,mu3,nu1,nu2):
	
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
		#print "coeff kronecker= ", kronecker( [l1,l2],[mu1,mu2,mu3],[nu1,nu2])
		pass
	
		#REGION P1 cap P2
		b1 = min(mu2,nu1/2)
		b0 = max(nu1-l2, (nu1-mu3)/2, nu1/3,mu3)		
		G1 = regionP1(l1,l2,mu1,mu2,mu3,nu1,nu2)
		G2 = regionP2(l1,l2,mu1,mu2,mu3,nu1,nu2)
		Gb = regionPbarra(l1,l2,mu1,mu2,mu3,nu1,nu2)
		GbT = regionPbarraT(l1,l2,mu1,mu2,mu3,nu1,nu2)

		show(G1,aspect_ratio=1.0)
		show(G2,aspect_ratio=1.0)
		
		#if (l1 >= mu1) and (mu1 >= l2) and (l2 >= mu2):
		#	print "interseccion:", floor(b1)-ceil(b0)+1
		#else:
		#	print "intersección: 0", "(fallan condiciones de orden)"
		
		show(Gb,aspect_ratio=1.0)
                print "Transformada de P-barra:"
		show(GbT,aspect_ratio=1.0)
                show(G1,aspect_ratio=1.0)
