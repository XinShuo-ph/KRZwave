
# coding: utf-8
from scipy.interpolate import interp1d
import numpy as np
from numpy import sin,cos,sqrt
#z1就是r，z2就是theta

def metric_KRZ(spin,   d1,   z1,   z2):
    r = z1;
    sqr = r*r;
    cuber = r*sqr;
    fourthr = sqr*sqr;

    th = z2;

    sinth = sin(th);
    costh = cos(th);
    sqsinth = sinth*sinth;
    sqcosth = costh*costh;

    sqspin = spin*spin;
    fourthspin = sqspin*sqspin;
    r0 = 1 + sqrt(1 - sqspin);
    sqr0 = r0*r0;
    cuber0 = sqr0*r0;
    fourthr0 = sqr0*sqr0;

    a20 = (2 * sqspin) / cuber0;
    a21 = -fourthspin / fourthr0 + 0;
    e0 = (2 - r0) / r0;
    k00 = sqspin / sqr0;
    k21 = fourthspin / fourthr0 - 2 * sqspin / cuber0 - 0;
    w00 = 2 * spin / sqr0;

    k22 = -sqspin / sqr0;
    k23 = sqspin / sqr0;

    N2 = (1 - r0 / r)*(1 - e0*r0 / r + (k00 - e0)*sqr0 / sqr + d1*cuber0 / cuber)+(a20*cuber0/cuber+a21*fourthr0/fourthr+k21*cuber0/cuber/(1+k22*(1-r0/r)/(1+k23*(1-r0/r))))*sqcosth;
    B = 1 + 0 * sqr0 / sqr + 0 * sqr0*sqcosth / sqr;
    Sigma = 1 + sqspin*sqcosth / sqr;
    W = (w00*sqr0 / sqr + 0 * cuber0 / cuber + 0 * cuber0*sqcosth / cuber) / Sigma;
    K2 = 1 + spin*W / r + (k00*sqr0 / sqr + k21*cuber0*sqcosth / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r)))) / Sigma;


    mn=np.zeros((4,4))
    mn[0][0] = -(N2 - W*W*sqsinth) / K2;
    mn[0][3] = -W*r*sqsinth;
    mn[1][1] = Sigma*B*B / N2;
    mn[2][2] = Sigma*sqr;
    mn[3][0] = mn[0][3];
    mn[3][3] = K2*sqr*sqsinth;

    return mn

def metric_KRZ_rderivatives_num(  spin,   d1,   z1,   z2):


    r = z1;
    theta = z2;
    dr = 0.001*r;


    mnm=metric_KRZ(spin, d1, r - dr, theta);
    mnp=metric_KRZ(spin, d1, r + dr, theta);
    
    rdmn=np.zeros((4,4))
    rdmn[0][0] = (mnp[0][0] - mnm[0][0])*0.5 / dr;
    rdmn[0][3] = (mnp[0][3] - mnm[0][3])*0.5 / dr;
    rdmn[1][1] = (mnp[1][1] - mnm[1][1])*0.5 / dr;
    rdmn[2][2] = (mnp[2][2] - mnm[2][2])*0.5 / dr;
    rdmn[3][0] = rdmn[0][3];
    rdmn[3][3] = (mnp[3][3] - mnm[3][3])*0.5 / dr;

    return rdmn
 

def metric_KRZ_thderivatives_num(  spin,   d1,   z1,   z2):


    r = z1;
    theta = z2;
    dtheta = 0.01;


    metric_KRZ(spin, d1, r, theta - dtheta, mnm);
    metric_KRZ(spin, d1, r, theta + dtheta, mnp);
    
    thdmn=np.zeros((4,4))
    thdmn[0][0] = (mnp[0][0] - mnm[0][0])*0.5 / dtheta;
    thdmn[0][3] = (mnp[0][3] - mnm[0][3])*0.5 / dtheta;
    thdmn[1][1] = (mnp[1][1] - mnm[1][1])*0.5 / dtheta;
    thdmn[2][2] = (mnp[2][2] - mnm[2][2])*0.5 / dtheta;
    thdmn[3][0] = thdmn[0][3];
    thdmn[3][3] = (mnp[3][3] - mnm[3][3])*0.5 / dtheta;

    return thdmn
 

def metric_KRZ_rderivatives(  spin,   d1,   z1,   z2) :

    r = z1;
    th = z2;

    sqr = r*r;
    cuber = r*sqr;
    fourthr = sqr*sqr;

    sinth = sin(th);
    costh = cos(th);
    sqsinth = sinth*sinth;
    sqcosth = costh*costh;

    sqspin = spin*spin;
    cubespin = sqspin*spin;
    fourthspin = sqspin*sqspin;

    r0 = 1 + sqrt(1 - sqspin);
    sqr0 = r0*r0;
    cuber0 = sqr0*r0;
    fourthr0 = sqr0*sqr0;

    a20 = (2 * sqspin) / cuber0;
    a21 = -fourthspin / fourthr0 + 0;
    e0 = (2 - r0) / r0;
    k00 = sqspin / sqr0;
    k21 = fourthspin / fourthr0 - 2 * sqspin / cuber0 - 0;
    w00 = 2 * spin / sqr0;

    k22 = -sqspin / sqr0;
    k23 = sqspin / sqr0;

    N2 = (1 - r0 / r)*(1 - e0*r0 / r + (k00 - e0)*sqr0 / sqr + d1*cuber0 / cuber) + (a20*cuber0 / cuber + a21*fourthr0 / fourthr + k21*cuber0 / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r))))*sqcosth;
    B = 1 + 0 * sqr0 / sqr + 0 * sqr0*sqcosth / sqr;
    Sigma = 1 + sqspin*sqcosth / sqr;
    W = (w00*sqr0 / sqr + 0 * cuber0 / cuber + 0 * cuber0*sqcosth / cuber) / Sigma;
    K2 = 1 + spin*W / r + (k00*sqr0 / sqr + k21*cuber0*sqcosth / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r)))) / Sigma;
    #2018-4-4
    #the oringinal N2 was wrong, see the correct KRZ metric

    rderN2 = (1 - r0 / r)*((2 - r0) / sqr - (2 * (-e0 + k00)*sqr0) / cuber - (3 * d1*cuber0) / fourthr) +    (r0*(1 - (2 - r0) / r + ((-e0 + k00)*sqr0) / sqr + (d1*cuber0) / cuber)) / sqr    +sqcosth*(-3*a20*cuber0/fourthr-4*a21*fourthr0/fourthr/r-3*k21*cuber0/fourthr/ (   1    +      k22*(1-r0/r)/(1+k23*(1-r0/r))  )-     k21*cuber0/cuber/  (1+  k22*(1-r0/r)/(1+k23*(1-r0/r ))  ) / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r))  )  *  ( k22*r0/sqr/(1+  k23*(1-r0/r)  ) - k22*(1-r0/r) / (1 + k23* (1 - r0 / r)) / (1 + k23* (1 - r0 / r)) *k23*r0/sqr  )    )
    #2018-4-4 modified with the correct metric expression

    '''+
    ((-3 * (-0 + fourthspin / fourthr0)*cuber0) / fourthr - (4 * a21*fourthr0) / np.power(r, 5))*sqcosth''';
    rderB = (-2 * 0 * sqr0) / cuber - (2 * 0 * sqr0*sqcosth) / cuber;
    rderSigma = (-2 * sqcosth*sqspin) / cuber;
    rderW = (2 * sqcosth*((0 * np.power(r0, 3)) / np.power(r, 3) + (2 * spin) / np.power(r, 2) + (0 * np.power(r0, 3)*sqcosth) / np.power(r, 3))*sqspin) /    (np.power(r, 3)*np.power(1 + (sqcosth*sqspin) / np.power(r, 2), 2)) +    ((-3 * 0 * np.power(r0, 3)) / np.power(r, 4) - (4 * spin) / np.power(r, 3) - (3 * 0 * np.power(r0, 3)*sqcosth) / np.power(r, 4)) /    (1 + (sqcosth*sqspin) / np.power(r, 2));
    '''	  rderK2 = (2 * cubespin*sqcosth*((0*cuber0) / cuber + (2 * spin) / sqr + (0*cuber0*sqcosth) / cuber)) /
    (fourthr*np.power(1 + (sqcosth*sqspin) / sqr, 2)) +
    (2 * sqcosth*sqspin*((k21*cuber0*sqcosth) / cuber + sqspin / sqr)) /
    (cuber*np.power(1 + (sqcosth*sqspin) / sqr, 2)) +
    (spin*((-3 *  0*cuber0) / fourthr - (4 * spin) / cuber - (3 * 0*cuber0*sqcosth) / fourthr)) /
    (r*(1 + (sqcosth*sqspin) / sqr)) - (spin*
    ((0*cuber0) / cuber + (2 * spin) / sqr + (0*cuber0*sqcosth) / cuber)) /
    (sqr*(1 + (sqcosth*sqspin) / sqr)) +
    ((-3 * k21*cuber0*sqcosth) / fourthr - (2 * sqspin) / cuber) / (1 + (sqcosth*sqspin) / sqr);'''
    rderK2 = spin / r*rderW - spin*W / sqr - rderSigma / Sigma / Sigma * (k00*sqr0 / sqr + k21*cuber0*sqcosth / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r))))    + (-2 * k00 * sqr0 / cuber - 3 * k21 * cuber0 * sqcosth / fourthr / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r)))    - k21 * sqcosth * cuber0 / cuber * (r0 * k22 / sqr / (1 + k23*(1 - r0 / r)) / (1 + k23*(1 - r0 / r))) / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r))) / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r)))) / Sigma;#2017-10-26

    rdmn=np.zeros((4,4))
    rdmn[0][0] = (-rderN2 + 2 * rderW*sqsinth*W) / K2 - (rderK2*(-N2 + sqsinth*np.power(W, 2))) / np.power(K2, 2);

    rdmn[1][1] = (np.power(B, 2)*rderSigma) / N2 + (2 * B*rderB*Sigma) / N2 - (np.power(B, 2)*rderN2*Sigma) / np.power(N2, 2);

    rdmn[2][2] = np.power(r, 2)*rderSigma + 2 * r*Sigma;

    rdmn[3][3] = 2 * K2*r*sqsinth + np.power(r, 2)*rderK2*sqsinth;

    rdmn[0][3] = -(r*rderW*sqsinth) - sqsinth*W;

    rdmn[3][0] = rdmn[0][3];

    return rdmn


def metric_KRZ_thderivatives(  spin,   d1,   z1,   z2) :
    r = z1;
    th = z2;

    sqr = r*r;
    cuber = r*sqr;
    fourthr = sqr*sqr;

    sinth = sin(th);
    costh = cos(th); #2018-8-19 之前算cosine居然是开根号。。。不知道以前的raytracing是不是也出过这种问题
    sqsinth = sinth*sinth;
    sqcosth = costh*costh;

    sqspin = spin*spin;
    cubespin = sqspin*spin;
    fourthspin = sqspin*sqspin;

    r0 = 1 + sqrt(1 - sqspin);
    sqr0 = r0*r0;
    cuber0 = sqr0*r0;
    fourthr0 = sqr0*sqr0;

    a20 = (2 * sqspin) / cuber0;
    a21 = -fourthspin / fourthr0 + 0;
    e0 = (2 - r0) / r0;
    k00 = sqspin / sqr0;
    k21 = fourthspin / fourthr0 - 2 * sqspin / cuber0 - 0;
    w00 = 2 * spin / sqr0;

    k22 = -sqspin / sqr0;#2017-10-26
    k23 = sqspin / sqr0;#2017-10-26

    N2 = (1 - r0 / r)*(1 - e0*r0 / r + (k00 - e0)*sqr0 / sqr + d1*cuber0 / cuber) + (a20*cuber0 / cuber + a21*fourthr0 / fourthr + k21*cuber0 / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r))))*sqcosth;
    B = 1 + 0 * sqr0 / sqr + 0 * sqr0*sqcosth / sqr;
    Sigma = 1 + sqspin*sqcosth / sqr;
    W = (w00*sqr0 / sqr + 0 * cuber0 / cuber + 0 * cuber0*sqcosth / cuber) / Sigma;
    K2 = 1 + spin*W / r + (k00*sqr0 / sqr + k21*cuber0*sqcosth / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r)))) / Sigma;
    #2018-4-4

    thderN2 = -2 * costh*sinth* (a20*cuber0 / cuber + a21*fourthr0 / fourthr + k21*cuber0 / cuber / (1 + k22* (1 - r0 / r) / (1 + k23 * (1 - r0 / r))));
    #2018-4-4 modified with the correct metric expression
    ''' -2 * costh*(((-0 + fourthspin / fourthr0)*cuber0) / cuber + (a21*fourthr0) / fourthr)*sinth ''';
    thderB = (-2 * costh * 0 * sqr0*sinth) / sqr;
    thderSigma = (-2 * costh*sinth*sqspin) / sqr;
    thderW = (2 * costh*sinth*((0 * cuber0) / cuber + (sqcosth * 0 * cuber0) / cuber + (2 * spin) / sqr)*    sqspin) / (sqr*np.power(1+(sqcosth*sqspin)/sqr, 2)) -    (2 * costh * 0 * cuber0*sinth) / (cuber*(1 + (sqcosth*sqspin) / sqr));
    '''  thderK2 = (2 * costh*cubespin*sinth*((0*cuber0) / cuber + (sqcosth*0*cuber0) / cuber +
    (2 * spin) / sqr)) / (cuber*np.power(1 + (sqcosth*sqspin) / sqr, 2)) +
    (2 * costh*sinth*sqspin*((sqcosth*k21*cuber0) / cuber + sqspin / sqr)) /
    (sqr*np.power(1 + (sqcosth*sqspin) / sqr, 2)) -
    (2 * costh*k21*cuber0*sinth) / (cuber*(1 + (sqcosth*sqspin) / sqr)) -
    (2 * costh*0*cuber0*sinth*spin) / (fourthr*(1 + (sqcosth*sqspin) / sqr));'''

    thderK2 = spin / r * thderW - thderSigma / Sigma / Sigma  * (k00*sqr0 / sqr + k21*cuber0*sqcosth / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r))))    + (-2* costh * sinth *  k21 * cuber0 / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r)))) / Sigma;#2017-10-26

    thdmn=np.zeros((4,4))
    
    thdmn[0][0] = (-thderN2 + 2 * np.power(sinth, 2)*thderW*W + 2 * costh*sinth*np.power(W, 2)) / K2 -    (thderK2*(-N2 + np.power(sinth, 2)*np.power(W, 2))) / np.power(K2, 2);
    thdmn[1][1] = (2 * B*Sigma*thderB) / N2 - (np.power(B, 2)*Sigma*thderN2) / np.power(N2, 2) + (np.power(B, 2)*thderSigma) / N2;
    thdmn[2][2] = np.power(r, 2)*thderSigma;
    thdmn[3][3] = 2 * costh*K2*np.power(r, 2)*sinth + np.power(r, 2)*np.power(sinth, 2)*thderK2;
    thdmn[0][3] = -(r*np.power(sinth, 2)*thderW) - 2 * costh*r*sinth*W;
    thdmn[3][0] = thdmn[0][3];


    return thdmn


def metric_KRZ_inverse(  spin,   d1,   z1,   z2) :
    r = z1;
    th = z2;

    sqr = r*r;
    cuber = r*sqr;
    fourthr = sqr*sqr;

    sinth = sin(th);
    costh = sqrt(1 - sinth*sinth);
    sqsinth = sinth*sinth;
    sqcosth = costh*costh;

    sqspin = spin*spin;
    cubespin = sqspin*spin;
    fourthspin = sqspin*sqspin;

    r0 = 1 + sqrt(1 - sqspin);
    sqr0 = r0*r0;
    cuber0 = sqr0*r0;
    fourthr0 = sqr0*sqr0;

    a20 = (2 * sqspin) / cuber0;
    a21 = -fourthspin / fourthr0 + 0;
    e0 = (2 - r0) / r0;
    k00 = sqspin / sqr0;
    k21 = fourthspin / fourthr0 - 2 * sqspin / cuber0 - 0;
    w00 = 2 * spin / sqr0;

    k22 = -sqspin / sqr0;#2017-10-26
    k23 = sqspin / sqr0;#2017-10-26
    '''
    N2 = (1 - r0 / r)*(1 - e0*r0 / r + (k00 - e0)*sqr0 / sqr + d1*cuber0 / cuber) + ((k21 + a20)*cuber0 / cuber + a21*fourthr0 / fourthr)*sqcosth;
    B = 1 + 0 * sqr0 / sqr + 0 * sqr0*sqcosth / sqr;
    Sigma = 1 + sqspin*sqcosth / sqr;
    W = (w00*sqr0 / sqr + 0 * cuber0 / cuber + 0 * cuber0*sqcosth / cuber) / Sigma;
    K2 = 1 + spin*W / r + (k00*sqr0 / sqr + k21*cuber0*sqcosth / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r)))) / Sigma;#2017-10-26
    '''

    N2 = (1 - r0 / r)*(1 - e0*r0 / r + (k00 - e0)*sqr0 / sqr + d1*cuber0 / cuber) + (a20*cuber0 / cuber + a21*fourthr0 / fourthr + k21*cuber0 / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r))))*sqcosth;
    B = 1 + 0 * sqr0 / sqr + 0 * sqr0*sqcosth / sqr;
    Sigma = 1 + sqspin*sqcosth / sqr;
    W = (w00*sqr0 / sqr + 0 * cuber0 / cuber + 0 * cuber0*sqcosth / cuber) / Sigma;
    K2 = 1 + spin*W / r + (k00*sqr0 / sqr + k21*cuber0*sqcosth / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r)))) / Sigma;#2017-10-26


    invg=np.zeros((4,4))
    invg[0][0] = -(K2 / N2);
    invg[1][1] = N2 / (np.power(B, 2)*Sigma);
    invg[2][2] = 1 / (sqr*Sigma);
    invg[3][3] = (-np.power(W, 2) + N2 / np.power(sin(th), 2)) / (K2*N2*sqr);
    invg[0][3] = -(W / (N2*r));
    invg[3][0] = invg[0][3];


    '''  g[4][4], gg;
    metric_KRZ(spin, d1, z1, z2, g);
    #metric_KRZ_inverse(spin, d1,  w1, w2, invg);
    gg = g[0][0] * g[3][3] - g[0][3] * g[0][3];
    if (std::fabs(gg) > 1e10) :
    invg[0][0] = g[3][3] / gg;
    invg[0][3] = -g[0][3] / gg;
    invg[1][1] = 1 / g[1][1];
    invg[2][2] = 1 / g[2][2];
    invg[3][0] = invg[0][3];
    invg[3][3] = g[0][0] / gg;
    '''
    return invg

def getwave(filename,THETA,PHI):
#读入文件名，和观测角，输出引力波
    try:
        index, tau,t,r,th,phi,ut,ur,uth,uphi,F_t,F_r,F_th,F_phi=np.loadtxt(filename,unpack=True)
    except:
        print(filename+'  does not exist')
        quit()

    #qseudo_flat spacetime
    x=[];
    y=[];
    z=[];
    t_tau_dot=[]
    z_tau_dot=[]
    y_tau_dot=[]
    x_tau_dot=[]
    z_t_dot=[]
    y_t_dot=[]
    x_t_dot=[]
    vr_tau_dot=[]
    vth_tau_dot=[]
    vphi_tau_dot=[]
    vx_tau_dot=[]
    vy_tau_dot=[]
    vz_tau_dot=[]
    x_t_2dot=[]
    y_t_2dot=[]
    z_t_2dot=[]


    for i in np.arange(index.size):
        x.append(r[i]*np.sin(th[i])*np.cos(phi[i]));
        y.append(r[i]*np.sin(th[i])*np.sin(phi[i]));
        z.append(r[i]*np.cos(th[i]));
        t_tau_dot.append(ut[i])
        x_tau_dot.append(ur[i]*np.sin(th[i])*np.cos(phi[i]) + r[i]*np.cos(th[i])*np.cos(phi[i])*uth[i] - r[i]*np.sin(th[i])*np.sin(phi[i])*uphi[i] )
        y_tau_dot.append(ur[i]*np.sin(th[i])*np.sin(phi[i]) + r[i]*np.cos(th[i])*np.sin(phi[i])*uth[i] + r[i]*np.sin(th[i])*np.cos(phi[i])*uphi[i] )
        z_tau_dot.append(ur[i]*np.cos(th[i]) - r[i]*np.sin(th[i])*uth[i])
        x_t_dot.append(x_tau_dot[i]/t_tau_dot[i])
        y_t_dot.append(y_tau_dot[i]/t_tau_dot[i])
        z_t_dot.append(z_tau_dot[i]/t_tau_dot[i])

        vr_tau_dot.append( (F_r[i]*t_tau_dot[i]-ur[i]*F_t[i])/t_tau_dot[i]/t_tau_dot[i] )
        vth_tau_dot.append( (F_th[i]*t_tau_dot[i]-uth[i]*F_t[i])/t_tau_dot[i]/t_tau_dot[i] )
        vphi_tau_dot.append( (F_phi[i]*t_tau_dot[i]-uphi[i]*F_t[i])/t_tau_dot[i]/t_tau_dot[i] )

        vx_tau_dot.append( vr_tau_dot[i]*np.sin(th[i])*np.cos(phi[i]) + ur[i]/ut[i]*np.cos(th[i])*np.cos(phi[i])*uth[i] - ur[i]/ut[i]*np.sin(th[i])*np.sin(phi[i])*uphi[i]\
             + ur[i]*cos(th[i])*cos(phi[i])*uth[i]/ut[i] - r[i]*sin(th[i])*cos(phi[i])*uth[i]/ut[i]*uth[i] -r[i]*cos(th[i])*sin(phi[i])*uth[i]/ut[i]*uphi[i] +r[i]*cos(th[i])*cos(phi[i])*vth_tau_dot[i]  \
             - ur[i]*sin(th[i])*sin(phi[i])*uphi[i]/ut[i] - r[i]*cos(th[i])*sin(phi[i])*uphi[i]/ut[i]*uth[i] - r[i]*sin(th[i])*cos(phi[i])*uphi[i]/ut[i]*uphi[i] - r[i]*sin(th[i])*sin(phi[i])*vphi_tau_dot[i])

        vy_tau_dot.append( vr_tau_dot[i]*np.sin(th[i])*np.sin(phi[i]) + ur[i]/ut[i]*np.cos(th[i])*np.sin(phi[i])*uth[i] + ur[i]/ut[i]*np.sin(th[i])*np.cos(phi[i])*uphi[i]\
             + ur[i]*cos(th[i])*sin(phi[i])*uth[i]/ut[i] - r[i]*sin(th[i])*sin(phi[i])*uth[i]/ut[i]*uth[i] +r[i]*cos(th[i])*cos(phi[i])*uth[i]/ut[i]*uphi[i] +r[i]*cos(th[i])*sin(phi[i])*vth_tau_dot[i]  \
             + ur[i]*sin(th[i])*cos(phi[i])*uphi[i]/ut[i] + r[i]*cos(th[i])*cos(phi[i])*uphi[i]/ut[i]*uth[i] - r[i]*sin(th[i])*sin(phi[i])*uphi[i]/ut[i]*uphi[i] + r[i]*sin(th[i])*cos(phi[i])*vphi_tau_dot[i])

        vz_tau_dot.append( vr_tau_dot[i]*cos(th[i]) -ur[i]/ut[i]*sin(th[i])*uth[i] \
                         -ur[i]*sin(th[i])*uth[i]/ut[i] -r[i]*cos(th[i])*uth[i]/ut[i]*uth[i] - r[i]*sin(th[i])*vth_tau_dot[i] )

        x_t_2dot.append(vx_tau_dot[i]/ut[i])
        y_t_2dot.append(vy_tau_dot[i]/ut[i])
        z_t_2dot.append(vz_tau_dot[i]/ut[i])

    #四极矩算法，在trace-reversed gauge的metric

    hbar_xx=[]
    hbar_yy=[]
    hbar_zz=[]
    hbar_xy=[]
    hbar_yz=[]
    hbar_xz=[]
    for i in np.arange(index.size):
        hbar_xx.append(4*(x_t_dot[i]*x_t_dot[i]+x[i]*x_t_2dot[i]))
        hbar_yy.append(4*(y_t_dot[i]*y_t_dot[i]+y[i]*y_t_2dot[i]))
        hbar_zz.append(4*(z_t_dot[i]*z_t_dot[i]+z[i]*z_t_2dot[i]))
        hbar_xy.append(2*(y[i]*x_t_2dot[i]+y_t_2dot[i]*x[i]+2*y_t_dot[i]*x_t_dot[i]))
        hbar_yz.append(2*(y[i]*z_t_2dot[i]+y_t_2dot[i]*z[i]+2*y_t_dot[i]*z_t_dot[i]))
        hbar_xz.append(2*(z[i]*x_t_2dot[i]+z_t_2dot[i]*x[i]+2*z_t_dot[i]*x_t_dot[i]))

    #由trace-reversed gauge转换到transverse traceless gauge

    hTT_TT=[]
    hTT_PP=[]
    hTT_TP=[]
    hTT_plus=[]
    hTT_cross=[]

    for i in np.arange(index.size):


        hTT_TT.append( np.cos(THETA)*np.cos(THETA)* (hbar_xx[i]*np.cos(PHI)*np.cos(PHI) + hbar_xy[i]*np.sin(2*PHI) + hbar_yy[i]*np.sin(PHI)*np.sin(PHI) )  +  hbar_zz[i]*np.sin(THETA)*np.sin(THETA)  -  np.sin(2*THETA)* (hbar_xz[i]*np.cos(PHI)+hbar_yz[i]*np.sin(PHI))  )
        hTT_TP.append( np.cos(THETA)* (-0.5*hbar_xx[i]*np.sin(2*PHI) + hbar_xy[i]*np.cos(2*PHI) + 0.5*hbar_yy[i]*np.sin(2*PHI))  +  np.sin(THETA)* (hbar_xz[i]*np.sin(PHI)-hbar_yz[i]*np.cos(PHI)) )
        hTT_PP.append( hbar_xx[i]*np.sin(PHI)*np.sin(PHI)  -  hbar_xy[i]*np.sin(2*PHI)  +  hbar_yy[i]*np.cos(PHI)*np.cos(PHI) )
        hTT_plus.append(0.5*(hTT_TT[i]-hTT_PP[i]))
        hTT_cross.append(hTT_TP[i])

    #注意上面算出来的h还要*mu（mass ratio）/R（观测距离，也以M为单位）才是真的strain
    #发现一个小问题，上面定义的数据类型大部分都是list，但是array才比较好用
    #还要注意一点几何单位制和SI单位的转换

    ########转换单位
    Grav=6.674e-11 #引力常数
    clight=2.998e8 #光速
    Msol=1.989e30  #太阳质量，以千克做单位

    M=1e6 # clight*clight*clight/Grav/Msol/1 #中心天体质量，以太阳质量为单位

    #把时间转换成秒
    t_sec=t*M*Msol*Grav/clight/clight/clight
    dt=t_sec[1]-t_sec[0]

    #把pc距离转换成M为单位
    R_pc=5e9  #以pc为单位的观测距离
    R=R_pc*3.0857e16*clight*clight/Grav/M/Msol  #以中心天体质量为单位的，长度米与中心天体质量的换算是 1m/kg = clight*clight/G

    #小天体的质量
    mu=1e-5 #应该是以中心天体质量为单位的

    hTT_plus_true=np.array(hTT_plus)*mu/R
    hTT_cross_true=np.array(hTT_cross)*mu/R

    ########用于计算的波形，plus作为实部，cross作为虚部

    return t,hTT_plus_true+hTT_cross_true*1j
    
def bracket(mydata,mytemp,dt,fnoise=[-1,-1],Snoise=[-1,-1]):
#算mydata和mytemp的内积
#fnoise-Snoise是噪声的功率谱
    
    #两个时间序列长度要一样，不一样就取最小的
    if len(mydata) != len(mytemp):
        print("inner product: length not match!")

    #采样频率
    fs=1/dt

    #做傅里叶变换，注意要除以采样频率才是真的amplitude
    mydata_fft=np.fft.fft(mydata)/len(mydata)
    mytemp_fft=np.fft.fft(mytemp)/len(mydata)
    #变换后的频率序列,注意，只有前一半是正的频率
    freq=np.fft.fftfreq(min( len(mydata),len(mytemp)),dt)
    
    #默认为LISA noise
    if fnoise[0]==-1:
        ##########LISA noise, reference: https://arxiv.org/abs/gr-qc/0607007v2
        u=2*np.pi*freq*50/3 #见reference（36）上面一段
        Sn=[]  #LISA noise
        for i in np.arange(freq.size/2):
            i=int(i)
            if i==0:
                Sn.append(1e10)
            elif(u[i]<0.25):
                Sn.append(8.08e-48/((2*np.pi*freq[i])**4) +5.52e-41 )
            else :
                Sn.append( (2.88e-48/((2*np.pi*freq[i])**4) +5.52e-41 ) *u[i]*u[i]/ ( (1+cos(u[i])*cos(u[i]) )*(1.0/3.0-2.0/u[i]/u[i]) + sin(u[i])**2 + 4*sin(u[i])*cos(u[i])/(u[i]**3) ) )

    else:
        try:#如果输入了噪声功率谱,先直接interpolate
            noisefit=interp1d(fnoise,Snoise)
            Sn=noisefit(freq[0:int(freq.size/2)])#如果输入的范围不够在这行会报错
            
        except:#如果报错说明输入的功率谱范围覆盖不了信号傅里叶变换后的功率谱，那么把输入功率谱的左右两侧用直线外推
            leftfit=np.poly1d(np.polyfit(fnoise[0:int(freq.size/8)],Snoise[0:int(freq.size/8)],1))
            rightfit=np.poly1d(np.polyfit(fnoise[int(freq.size/2)-int(freq.size/8):int(freq.size/2)],Snoise[int(freq.size/2)-int(freq.size/8):int(freq.size/2)],1))
            Sn=[]
            for i in np.arange(int(freq.size/2)):
                if freq[i]<fnoise[0]:
                    Sn.append(leftfit(freq[i]))
                elif freq[i]>fnoise[-1]:
                    Sn.append(rightfit(freq[i]))
                else:
                    Sn.append(noisefit(freq[i]))
                    
    ###########SNR
    SNR=0
    df=freq[1]-freq[0]
    for i in np.arange(int(freq.size/2)):
        i=int(i)
        SNR=SNR+(( mydata_fft[i]*np.conjugate(mytemp_fft[i]) + np.conjugate(mydata_fft[i])*mytemp_fft[i])/(Sn[i]))*df
    
    return np.abs(SNR)

def getfreq_fromtraj(tau,r,phi):
    #由序列获得orbital frequency
    p=np.mean(r)
    omgr=[]
    omgphi=[]
    indr=[]
    phi=phi-phi[0]
    n=1
    for i in np.arange(tau.size-1):
        if(r[i]>p and r[i+1]<=p):
            indr.append(i)

    for ii in np.arange(len(indr)):
        if ii==0:
            continue
        omgr.append(2*np.pi/(tau[indr[ii]]-tau[indr[ii-1]]))
        omgphi.append((phi[indr[ii]]-phi[indr[ii-1]])/(tau[indr[ii]]-tau[indr[ii-1]]))

    omgr=np.array(omgr)
    avgomgr=np.mean(omgr)
    avgomgphi=np.mean(np.array(omgphi))

    return avgomgr,avgomgphi
def getfreq_frommaxi(tau,r,phi):
    #由序列获得orbital frequency,这个序列必须是从r最大值开始的
    p=np.mean(r)
    omgr=[]
    omgphi=[]
    indr=[]
    phi=phi-phi[0]
    n=1
    for i in np.arange(tau.size-1):
        if i==0:
            indr.append(i)
        elif(r[i]>r[i-1] and r[i+1]<r[i]):
            indr.append(i)

    for ii in np.arange(len(indr)):
        if ii==0:
            continue
        omgr.append(2*np.pi/(tau[indr[ii]]-tau[indr[ii-1]]))
        omgphi.append((phi[indr[ii]]-phi[indr[ii-1]])/(tau[indr[ii]]-tau[indr[ii-1]]))

    omgr=np.array(omgr)
    avgomgr=np.mean(omgr)
    avgomgphi=np.mean(np.array(omgphi))

    return avgomgr,avgomgphi
def getfreq_fromepa(e,p,spin):

    rmax=p/(1-e)
    rmin=p/(1+e)
    print('called')
    invgmin=metric_KRZ_inverse(spin,0,rmin,np.pi/2)
    invgmax=metric_KRZ_inverse(spin,0,rmax,np.pi/2)

    EoverL = ((invgmax[3][0] - invgmin[3][0]) + sqrt((invgmax[3][0] - invgmin[3][0]) *(invgmax[3][0] - invgmin[3][0]) - (invgmax[0][0] - invgmin[0][0])*(invgmax[3][3] - invgmin[3][3]))) / ( invgmax[0][0]-invgmin[0][0] );
    Lz = sqrt( (invgmax[3][0]-invgmin[3][0]) / ( EoverL*EoverL*( invgmin[3][0]*invgmax[0][0] - invgmax[3][0]*invgmin[0][0] )+ ( invgmin[3][0]*invgmax[3][3]- invgmax[3][0]*invgmin[3][3] )  )   );
    E=Lz*EoverL
    Q=0

    dchi=1e-8#积分的精度
    chi=np.linspace(dchi/2,np.pi-dchi/2,int(np.pi/dchi))

    #reference: https://arxiv.org/pdf/gr-qc/0202090.pdf (39)-(42)
    myJ=(1-E**2)*(1-e**2)+2*(1-E**2-(1-e**2)/p)*(1+e*np.cos(chi))+((1-E**2)*(3+e**2)/(1-e**2) -4/p +(spin**2*(1-E**2)+Lz**2+Q)*(1-e**2)/p**2 )*np.power(1+e*np.cos(chi),2)
    myH=1-2/p*(1+e*cos(chi))+(spin**2)/(p**2)*np.power(1+e*np.cos(chi),2)
    myG=Lz-2*(Lz-spin*E)/p*(1+e*np.cos(chi))

    myY=p**2*dchi*np.sum(np.multiply( np.power(1+e*np.cos(chi),-2) ,np.power(myJ,-0.5) ))
    myZ=dchi*np.sum(np.multiply(myG, np.multiply( np.power(myH,-1),np.power(myJ,-0.5) ) ))

    myomgr=np.pi*p/(1-e**2)/myY
    myomgphi=myZ/myY

    return myomgr,myomgphi