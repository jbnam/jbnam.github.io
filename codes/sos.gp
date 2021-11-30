default(parisize, "1G");

\p 76;

JB_zeta(s) = lfun(lfuncreate(1), s);
JB_zeta_der(s, m=1) = lfun(lfuncreate(1), s, m);

JB_lchi4(s) = lfun(lfuncreate(-4), s);
JB_lchi4_der(s, m=1) = lfun(lfuncreate(-4), s, m);

/* All these constants are defined on Section 9. List of Constants */
K = 0.764223653589220662990698731250092328116790541393409514721687;
c_Shanks = 0.5819486593172907979281498845023675593048328730717734215249;
JB_gamma = 0.5772156649015328606065120900824024310421593359399235988058;
JB_gamma1 = -0.07281584548367672486058637587490131913773633833433795259901;
JB_gamma2 = -0.00969036319287231848453038603521252935906580610134074988070;
JB_alpha1 = 0.245609584777314172388816626179062518433533782954926619588;
JB_alpha2 = -0.19625933903710052458751632882525623949324116970452959148357;
JB_alpha3 = 0.120808608458113395057401241932844490080511031943853683520;
JB_beta1 = 0.4574727063880649152906961328127586186195244768216850516;
JB_beta2 = 0.8296140548391727790649240430639518742186190145728626447;
JB_beta3 = 0.9829527142934445912581645562812494944560171298276995153;
JB_omega = log(2/(Pi^2)) + JB_alpha1 + JB_beta1;
JB_omega1 = (Pi^2)/12 + 3*((log(2))^2) - JB_gamma^2 - 2*JB_gamma1 + \
			(JB_omega^2)/4 + (JB_alpha1^2)/2 - JB_alpha2/2 + JB_beta2;
JB_omega2 = 2*JB_zeta(3) -2*JB_gamma^3 -6*JB_gamma*JB_gamma1 -3*JB_gamma2 \
			-6*(log(2)^3) -JB_alpha3/2 +3*JB_alpha1*JB_alpha2/2 -JB_alpha1^3 \
			-4*JB_beta3 +((log(2) +JB_alpha1 +JB_beta1)^3)/8 \
			+(log(2) +JB_alpha1 +JB_beta1)/2 +(log(2) + JB_alpha1 + JB_beta1)*((JB_alpha1^2)/2 \
			+ log(2)^2 -JB_alpha2/2 +JB_beta2)/2 -3*JB_omega*JB_omega1 + (JB_omega^3)/4;

/*lcqchi = read("JB_gp/D_5_chi.gp");*/

/* n mod 4 */
JB_mod4(n) = {return(bitand(n,3))};

/* alpha(x) = 1-K/sqrt(log(x)) */
JB_alpha(x) = return(1.-K/(log(x)^(0.5)));

/* H = -1/log(alpha(x)) */
JB_H(x) = return(-1./log(JB_alpha(x)));

/* initialise the primes = 3 mod 4 up to Max_p */
Max_p = 10^5;
lp = List();
forprime(p=3,Max_p,if(p%4 == 3, listput(lp,p)));

/* prod_{p=3mod4}(1-1/p^{s*2^j}) */
JB_ep_pow_2s(s,j)=
{
	r = 1.;
	t = s*(2.^j);
	for(n=1,length(lp), pp=lp[n]^t; r*=((pp-1.)/pp););
	return(r);
}

/*  prod_{p=3mod4}(1-1/p^{2s}) */
JB_ep_2s(s, num_j=4)=
{
	r = 1.;
	for(j=1,num_j,t=s*(2.^j);rr=JB_lchi4(t)/JB_zeta(t)/(1.-(2.^(-t)));r*=(rr^(1/(2^j))););
	return(r*(JB_ep_pow_2s(s, num_j+1)^(1/(2^num_j))));
}

/* 2sum_{p=3mod4}log(p)/(p^{s*2^{j}}-1) */
JB_ep_log_der_pow_2s(s,j) = 
{
	r = 0.;
	t = s*(2.^j);
	for(n=1,length(lp), pp=lp[n]^t; r+=(log(lp[n])/(pp-1.)));
	return(2.*r);
}

/* 2sum_{p=3mod4}log(p)/(p^(2s)-1) = logarithmic derivative of prod_{p=3mod4}(1-p^(-2s)) 
with respect to s using j times iterations */
JB_ep_2s_log_der(s,num_j=4)=
{
	r = 0.;
	for(j=1,num_j,t=s*(2.^j); r+=(JB_lchi4_der(t)/JB_lchi4(t)-JB_zeta_der(t)/JB_zeta(t)
		-log(2.)/((2.^t)-1.)););
	return(r+JB_ep_log_der_pow_2s(s,num_j+1));
}

/* -4sum_{p=3mod4}(log(p))^2p^(2s)/(p^(2s)-1)^2 = derivative of the logarithmic derivative of 
prod_{p=3mod4}(1-p^(-2s)) with respect to s using j times iterations */
JB_ep_2s_der_log_der(s,num_j=4)=
{
	r = 0.;
	for(j=1,num_j,t=2^j;st=s*t;zd=JB_zeta(st);zd1=JB_zeta_der(st,1);zd2=JB_zeta_der(st,2);
		ld=JB_lchi4(st);ld1=JB_lchi4_der(st,1);ld2=JB_lchi4_der(st,2);ef2=2.^st;
		r+=(t*(ld2/ld-(ld1/ld)^2-zd2/zd+(zd1/zd)^2+(log(2)^2)*ef2/((ef2-1.)^2))));
	r1 = 0.;
	t1 = 2.^(num_j+1);
	for(j=1,length(lp),p=lp[j]^t1;r1+=((log(lp[j])^2)*p/((p-1.)^2)););
	return(r-2.*r1*t1);
}

/* 8sum_{p=3mod4}(log(p))^3p^(2s)(p^(2s)+1)/(p^(2s)-1)^3 = 2nd derivative of the logarithmic derivative of 
prod_{p=3mod4}(1-p^(-2s)) with respect to s using j times iterations */
JB_ep_2s_2nd_der_log_der(s,num_j=4)=
{
	r = 0.;
	for(j=1,num_j,t=2^j;st=s*t;
		zd=JB_zeta(st);zd1=JB_zeta_der(st,1);zd2=JB_zeta_der(st,2);zd3=JB_zeta_der(st,3);
		ld=JB_lchi4(st);ld1=JB_lchi4_der(st,1);ld2=JB_lchi4_der(st,2);ld3=JB_lchi4_der(st,3);ef2=2.^st;
		r+=((t^2)*(ld3/ld-3.*ld1*ld2/(ld^2)+2.*((ld1/ld)^3)
			-zd3/zd+3.*zd1*zd2/(zd^2)-2.*((zd1/zd)^3)-(log(2)^3)*ef2*(ef2+1.)/((ef2-1.)^3))));
	r1 = 0.;
	t1 = 2.^(num_j+1);
	for(j=1,length(lp),p=lp[j]^t1;r1+=((log(lp[j])^3)*p*(p+1.)/((p-1.)^3)););
	return(r-2.*r1*(t1^2));
}

/* M(s) = (1-2^(-s)+2^(-2s))(L(s+1, chi_4)(1-2^(-(s+1)))prod_{p=3mod4}(1-p^(-2(s+1))))^(-1/2) */
JB_M(s) = 
{
	t = 1./(2.^s);
	return((1.-t+t^2)/((JB_lchi4(s+1.)*(1.-(t/2.))*JB_ep_2s(s+1.))^(1/2)));
}

/* logarithmic derivative of M(s) */
JB_M_log_der(s) = 
{
	t = 2.^s;
	return(log(2)*(t-2.)/(t^2-t+1.) -(JB_lchi4_der(s+1.)/JB_lchi4(s+1.) + log(2.)/(2.*t-1.) + JB_ep_2s_log_der(s+1.))/2.);
}

/* F(s) = zeta(s-1)M(s-1)sqrt((s-1)zeta(s))Gamma(s) */
JB_F(s) = return(JB_zeta(s-1.)*JB_M(s-1.)*(((s-1.)*JB_zeta(s))^(1/2))*gamma(s));

/* F_0(s,q) = A_0(s,q)F(s) = ((1-q^(-(s-1)))/(s-1))F(s) 
JB_F_0(s,q) = return((1.-q^(-(s-1.)))*JB_zeta(s-1.)*JB_M(s-1.)*(((s-1.)*JB_zeta(s))^(1/2))*gamma(s-1.));
*/

/* A_0(s, q) = (1-q^(-(s-1)))/(s-1) */ 
JB_A_0(s, q) = return((1.-q^(-(s-1.)))/(s-1.));

/* F_0(s,q) = F(s)A_0(s,q) = F(s)(1-q^(-(s-1)))/(s-1) */ 
JB_F_0(s, q) = return(JB_F(s)*JB_A_0(s, q));

/* logarithmic derivative of F(s) */
JB_F_log_der(s) = return(JB_zeta_der(s-1.)/JB_zeta(s-1.) + JB_M_log_der(s-1.) + (1./(s-1.) + JB_zeta_der(s)/JB_zeta(s))/2. + psi(s));

/* the Lucile's integrands in the integral form for S(q,v,H) */
JB_SqvH_integrand_Lucile(s, q, v, h) = 
{	
	if(v%q==0,return(JB_F(s)*((h/q)^(s-1.))*(JB_F_log_der(s)+log(h/q))/(abs(s-1.)^(1/2)));,return((1.-q^(s-1.))*JB_F(s)*((h/q)^(s-1.))/(abs(s-1.)^(3/2))););
}

/* the Chantal's integrands in the integral form for S(q,v,H) */
JB_SqvH_integrand_Chantal(s, q, v, h) = 
{	
	if(v%q==0,return(JB_F(s)*(h^(s-1.))*(JB_F_log_der(s)+log(h)-JB_A_0(s,q)/2.)/(abs(s-1.)^(1/2)));,return(JB_F_0(s,q)*(h^(s-1.))/(abs(s-1.)^(1/2))););
}

/* the Lucile's integrals in the integral form for S(q,v,H) */
JB_SqvH_integral_Lucile(q, v, h, prec=20) = 
{
	r = 1./Pi/(K^2);
	if(v%q!=0,r*=(1./2./eulerphi(q)));
	return(r*intnum(t=0.5+(1./10^prec), 1.-(1./10^prec),JB_SqvH_integrand_Lucile(t,q,v,h)));
}

/* the Chantal's integrals in the integral form for S(q,v,H) */
JB_SqvH_integral_Chantal(q, v, h, prec=20) = 
{
	r = 1./Pi/(K^2);
	if(v%q!=0,r*=(1./2./eulerphi(q)));
	return(r*intnum(t=0.5+(1./10^prec), 1.-(1./10^prec),JB_SqvH_integrand_Chantal(t,q,v,h)));
}

/* evaluation of a Dirichlet character chi^(2^j) at n */
JB_chi(n, vchi, j=0) =  return(exp(2.*Pi*I*chareval(vchi[1], vchi[2]*(2^j), n)));

/* List of pairs of (L(s,chiq), L(s,chiq4)) where chiq and chiqchi4 are Dirichlet characters
mod q and mod 4q, repectively, where q is a prime and congruent to 1 mod 4 such that chiq is
odd and chi4 is the primitive character mod 4 */
JB_G(q) = 
{
	lchiqchiq4 = List();
	q4 = 4*q;
	gq4 = znstar(q4,1);
	gq = znstar(q,1);
	for(j=1,q4,if(gcd(j,q4)==1,chiqchiq4=znconreychar(gq4,j);
		if(#znconreyconductor(gq4,chiqchiq4)==1,chiq=[chiqchiq4[1]]; 
			if(zncharisodd(gq,chiq)==1,listput(lchiqchiq4,[[gq,chiq],[gq4,chiqchiq4]]);););););
	return(lchiqchiq4);
}

/* L(0,chiq)L(1,chiq)^(1/2)L(1,chiqchiq4)^(-1/2)(1-chiq(2)-chiq(4))
(1-chiq(2)/2)^(-1/2)prod_{p=3mod4}(1-chi(p)^2p^(-2))^(-1/2) */
JB_Cqchi(vchiqchiq4, num_j=4)=
{
	r = lfun(lfuncreate(vchiqchiq4[1]),0)*(lfun(lfuncreate(vchiqchiq4[1]),1)^(1/2))
	*(1.-JB_chi(2,vchiqchiq4[1])+JB_chi(4,vchiqchiq4[1]))/((1.-JB_chi(2,vchiqchiq4[1])/2.)^(1/2))
	/(lfun(lfuncreate(vchiqchiq4[2]),1)^(1/2));
	r1 = 1.;
	r2 = 1.;
	for(j=1,num_j,t=2^j;if(t!=charorder(vchiqchiq4[1][1],vchiqchiq4[1][2]),
		tvq=[vchiqchiq4[1][1],vchiqchiq4[1][2]*t];
		tvqq4=[vchiqchiq4[2][1],[vchiqchiq4[2][2][1]*t,vchiqchiq4[2][2][2]]];
		r1*=((lfun(lfuncreate(tvqq4),t)/lfun(lfuncreate(tvq),t)/(1.-JB_chi(2,tvq)/(2.^t)))^(1/t));
		,r1*=JB_ep_2s(t/2)^(2/t));return(r/(r1^(1/2))););
	t1 = 2^(num_j+1);
	for(j=1,length(lp),pp=lp[j]^t1;r2*=((pp-JB_chi(lp[j],[vchiqchiq4[1][1],vchiqchiq4[1][2]*t1]))/pp););
	r2 = r2^((1/2)^num_j);
	return(r/((r1*r2)^(1/2)));
}

/* G(s) in the integral in Theorem 2.4 */
JB_integrand_2_4(s, num_iter)=
	return(((abs(real(JB_zeta(s)))*JB_lchi4(s)/JB_ep_2s(s, num_iter)/(1.-(0.5)^s))^(0.5))/s);

/* The integral in Theorem 2.4 */
JB_integral_2_4(x, num_iterations=4)=
	return(intnum(t=0.5+1.E-10, 1., (x^t)*JB_integrand_2_4(t, num_iterations))/Pi);

/* f(v;q) = -1/2 if v = 0 and (q-2v)/2/q = 1/2 -v/q otherwise. */
JB_f(q,v) = {vv=v%q; if(vv==0, vv=q); return(1./2-vv/q)};

/* C_{a,b}: C_{q\chi} are in l */
JB_Cqab(q, v) = 
{
	vv=v%q;
	if(vv==0,return(0.));
	y = 0.;
	for(j=2,length(lcqchi),y+=lcqchi[j][2]/lcqchi[j][1][vv+1];);
	return(real(y)/2./eulerphi(q)/(K^2));
}
