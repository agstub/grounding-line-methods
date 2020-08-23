# These functions produce sample plots for the test problem described in the README.
import sys
sys.path.insert(0, './source')

import matplotlib.pyplot as plt
import numpy as np
from geometry import interface, bed
from params import tol,X_fine,nt_per_year,t_final,Lngth,Hght,sea_level
from hydrology import *

#-------------------------------------------------------------------------------
#----------------------Plotting for test problem--------------------------------
def testplot_2(Gamma_s,Gamma_h,x_left,x_right):
    X = X_fine

    B = np.array([25,85,205,249])

    plt.figure(figsize=(12,8))
    j=1

    labels = [r'(a) lowstand',r'(b) filling',r'(c) highstand',r'(d) draining']

    for i in B:
        str1 = '22'+str(j)
        plt.subplot(int(str1))


        plt.annotate(labels[j-1],xy=(0.1,11.5),fontsize=20)
        plt.plot(X/1000,Gamma_h[:,i]-Hght+ 10,color='royalblue',linewidth=1.5,label=r'$h$')

        plt.fill_between(X/1000,y1=Gamma_s[:,i], y2=Gamma_h[:,i]-Hght+ 10,facecolor='aliceblue',alpha=1.0)
        plt.fill_between(X/1000,bed(X),Gamma_s[:,i],facecolor='slateblue',alpha=0.5)
        plt.fill_between(X/1000,-13*np.ones(np.size(X)),bed(X),facecolor='burlywood',alpha=1.0)

        plt.plot(X/1000,bed(X),color='k',linewidth=1.5,label=r'$\beta$')
        plt.plot(X[Gamma_s[:,i]-bed(X)>tol]/1000,Gamma_s[:,i][Gamma_s[:,i]-bed(X)>tol],color='crimson',linewidth=1.5,label=r'$s>\beta$')

        plt.plot(np.array([x_left[i]/1000]),np.array(np.min(bed(X))+1-2),marker='^',color='k',linestyle='None',markersize=10,label=r'$x_\pm$')
        plt.plot(np.array([x_right[i]/1000]),np.array(np.min(bed(X))+1-2),marker='^',markersize=10,color='k')

        if j == 1:
            plt.plot(X/1000,Gamma_h[:,0]-Hght+ 10,color='k',linestyle='--',linewidth=1.5)
            plt.plot(X[Gamma_s[:,0]-bed(X)>tol]/1000,Gamma_s[:,0][Gamma_s[:,0]-bed(X)>tol],color='k',linestyle='--',linewidth=1.5)



        if j == 1 or j==3:
            plt.yticks([-10,-5,0,7.5,10,12.5],[r'$-10$',r'$-5$',r'$0$',r'$h_0-2.5$',r'$h_0$',r'$h_0+2.5$'],fontsize=16)
            plt.ylabel(r'$z$ (m)',fontsize=20)
        else:
            plt.yticks([-10,-5,0,7.5,10,12.5],[r'',r'',r'',r'',r'',r''],fontsize=16)
            plt.gca().get_yaxis().set_ticklabels([])

        if j >2:
            plt.xlabel(r'$x$ (km)',fontsize=20)
            plt.xticks([0,1,2,3,4,5,6,7,8,9,10],['-5','-4','-3','-2','-1','0','1','2','3','4','5'],fontsize=16)
        else:
            plt.xticks([0,1,2,3,4,5,6,7,8,9,10],['','','','','','','','','','',''],fontsize=16)

        plt.yticks(fontsize=16)
        plt.ylim(np.min(bed(X))-2,14.0,8)
        plt.xlim(0,Lngth/1000)

        j+=1

    plt.tight_layout()
    lgd = plt.legend(fontsize=20,bbox_to_anchor=(0.55, -0.225),ncol=4)
    plt.savefig('laketest2', bbox_extra_artists=(lgd,),bbox_inches='tight')
    plt.close()

#-------------------------------------------------------------------------------

def testplot_1(s_mean,h_mean,lake_vol,x_right,x_left,dPw):

    t = np.linspace(0,t_final,num=int((t_final/3.154e7)*nt_per_year))

    B = np.array([25,85,220,249])

    plt.figure(figsize=(13,6))
    plt.subplot(221)
    plt.plot(t/3.154e7,Vol(t,1),'k',linewidth=3,label=r'true')
    plt.plot(t/3.154e7,lake_vol/lake_vol[0],'seagreen',linewidth=3,linestyle='--',label=r'computed')
    plt.plot(t[B]/3.154e7,Vol(t[B],1),'kD',markersize=10)
    plt.ylabel(r'$\mathcal{V}\,/\,\mathcal{V}_0$',fontsize=20)
    plt.gca().get_xaxis().set_ticklabels([])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(fontsize=16,loc='lower right',framealpha=1.0)


    plt.subplot(222)
    p1, = plt.plot(t/3.154e7,s_mean-s_mean[0],'crimson',label=r'$\Delta\bar{s}$',linewidth=3)
    p2, = plt.plot(t/3.154e7,h_mean-h_mean[0],'royalblue',label=r'$\Delta\bar{h}$',linewidth=3)
    plt.ylabel(r'$z$ (m)',fontsize=20)
    plt.gca().get_xaxis().set_ticklabels([])
    plt.gca().yaxis.tick_right()
    plt.gca().yaxis.set_label_position("right")
    plt.ylim(-2,2)
    plt.xticks(fontsize=16)
    plt.yticks(np.arange(-2,2.1,1),fontsize=16)
    plt.legend(fontsize=16,loc='lower right',framealpha=1.0)

    plt.subplot(223)
    plt.fill_between(t/(3.154e7),y1=x_left/1000.0 -5, y2=x_right/1000.0 - 5,facecolor='mediumslateblue',alpha=0.5,label='lake water')
    plt.plot(t/3.154e7,x_right/1000.0 - 5,'o',markersize=3,color='k')
    plt.plot(t/3.154e7,x_left/1000.0 -5 ,'o',markersize=3,color='k')
    plt.gca().ticklabel_format(axis='x',style='sci')
    plt.gca().yaxis.major.formatter._useMathText = True
    plt.xlabel(r'$t$ (yr)',fontsize=20)
    plt.ylabel(r'$x_\pm$ (km)',fontsize=20)
    plt.xticks(fontsize=16)
    plt.yticks(np.arange(-3.0,3.5,1.0),fontsize=16)
    plt.annotate(r'lake water',xy=(1.25,-0.25),fontsize=20)
    plt.ylim(-4,4)

    plt.subplot(224)
    plt.axhline(y=0.0,linewidth=2,linestyle='-',color='navy')
    plt.plot(t/3.154e7,dPw,'ko',markersize=3)
    plt.gca().ticklabel_format(axis='x',style='sci')
    plt.gca().yaxis.major.formatter._useMathText = True
    plt.xlabel(r'$t$ (yr)',fontsize=20)
    plt.gca().yaxis.tick_right()
    plt.gca().yaxis.set_label_position("right")

    plt.ylabel(r'$\bar{p}_\mathrm{w}-p_\mathrm{o}$ (kPa)',fontsize=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.tight_layout()
    plt.savefig('laketest1', bbox_inches='tight')

#-------------------------------------------------------------------------------




#-----------------Plotting for problems in paper (Figs 2-5)---------------------
#-------------------------------------------------------------------------------

def paperplot_fig2(s_mean,h_mean,x_right,x_left,Gamma_h,Gamma_s):

    x_left = x_left/1000.0 - 0.5*Lngth/1000.0
    x_right = x_right/1000.0 - 0.5*Lngth/1000.0

    h_out = Gamma_h[-1,:]
    s_out = Gamma_s[-1,:]


    t = np.linspace(0,t_final,num=int((t_final/3.154e7)*nt_per_year))

    h_f = (1.0/0.917)*(sea_level+sl_change(t))+(1-1.0/0.917)*s_out

    B = 5*np.array([1495,1528,1563,1597]) #Indices corresponding to free surface plots

    plt.figure(figsize=(8,8))
    plt.subplot(311)
    plt.title(r'(a)',fontsize=20,loc='left')
    plt.plot(t/(3.154e7/12.0/30.0),s_mean-s_mean[0],'crimson',linewidth=3,linestyle='-',label=r'$\Delta \bar{s}$')
    plt.plot(t/(3.154e7/12.0/30.0),h_mean-h_mean[0],'royalblue',linewidth=3,linestyle='-',label=r'$\Delta \bar{h}$')
    plt.plot(t/(3.154e7/12.0/30.0),sl_change(t),'k',linewidth=3,linestyle='-',label=r'$\Delta\ell$')
    plt.plot(t[B]/(3.154e7/12.0/30.0),sl_change(t)[B],'kD',markersize=10)
    plt.ylabel(r'$z$ (m)',fontsize=20)
    plt.gca().get_xaxis().set_ticklabels([])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(fontsize=16,bbox_to_anchor=(1.0,0.75))
    plt.xlim(5,7)

    plt.subplot(312)
    plt.title(r'(b)',fontsize=20,loc='left')
    plt.plot(t/(3.154e7/12.0/30.0),h_f-Hght,'k',linewidth=3,linestyle='-',label=r'$\Delta h_{\mathrm{out}}$')
    plt.plot(t/(3.154e7/12.0/30.0),h_out-Hght,'seagreen',linewidth=3,linestyle='--',label=r'$\Delta h_f$')
    plt.ylabel(r'$z$ (m)',fontsize=20)
    plt.gca().get_xaxis().set_ticklabels([])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(fontsize=16,bbox_to_anchor=(1.0,0.75))
    plt.xlim(5,7)



    plt.subplot(313)
    plt.title(r'(c)',fontsize=20,loc='left')
    plt.fill_between(t/(3.154e7/12.0/30.0),y1=x_left, y2=x_right,facecolor='rosybrown',alpha=1.0,label='partial \ncontact')
    plt.plot(t/(3.154e7/12.0/30.0),x_left,'o',color='k',markersize=5)
    plt.plot(t/(3.154e7/12.0/30.0),x_right,'o',color='k',markersize=5)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.ylabel(r'$x_\pm$ (km)',fontsize=20)
    plt.xlabel(r'$t$ (d)',fontsize=20)
    plt.legend(fontsize=16,markerscale=1.25,bbox_to_anchor=(1.0,0.25))
    plt.ylim(-4,4)
    plt.xlim(5,7)
    plt.tight_layout()
    plt.savefig('fig2', bbox_inches='tight')
    plt.close()
#-------------------------------------------------------------------------------

def paperplot_fig3(Gamma_s,Gamma_h,x_left,x_right):
    X = X_fine

    # Sample time indices
    B = 5*np.array([1495,1528,1563,1597])

    t = np.linspace(0,t_final,num=int((t_final/3.154e7)*nt_per_year))


    plt.figure(figsize=(12,8))
    j=1

    labels = [r'(a) low tide',r'(b) rising sea level',r'(c) high tide',r'(d) falling sea level']

    for i in B:
        str1 = '22'+str(j)
        plt.subplot(int(str1))

        plt.annotate(labels[j-1],xy=(0.1,17.5),fontsize=20)

        plt.plot(X/1000,Gamma_h[:,i]-Hght+ 10,color='royalblue',linewidth=1.5,label=r'$h$')

        plt.fill_between(X/1000,y1=Gamma_s[:,i], y2=Gamma_h[:,i]-Hght+ 10,facecolor='aliceblue',alpha=1.0)
        plt.fill_between(X/1000,bed(X),Gamma_s[:,i],facecolor='slateblue',alpha=0.5)
        plt.fill_between(X/1000,-13*np.ones(np.size(X)),bed(X),facecolor='burlywood',alpha=1.0)

        plt.plot(X/1000,bed(X),color='k',linewidth=1.5,label=r'$\beta$')

        plt.plot(X[(Gamma_s[:,i]-bed(X)>tol)&(X/1000.0>=x_right[i]/1000.0)]/1000,Gamma_s[:,i][(Gamma_s[:,i]-bed(X)>tol)&(X/1000.0>=x_right[i]/1000.0)],color='crimson',linewidth=1,label=r'$s>\beta$')
        plt.plot(X[(Gamma_s[:,i]-bed(X)>tol)&(X/1000.0<x_right[i]/1000.0)]/1000,Gamma_s[:,i][(Gamma_s[:,i]-bed(X)>tol)&(X/1000.0<x_right[i]/1000.0)],'o',color='crimson',markersize=1.5)

        plt.plot(np.array([x_left[i]/1000]),np.array(np.min(bed(X))+1-2),marker='^',color='k',linestyle='None',markersize=10,label=r'$x_\pm$')
        plt.plot(np.array([x_right[i]/1000]),np.array(np.min(bed(X))+1-2),marker='^',markersize=10,color='k')

        plt.axhline(y=0.02*sea_level+sl_change(t[i]),xmin=0.8,linestyle='--',color='seagreen',linewidth=1.5,label=r'$\ell$')
        plt.axhline(y=(1.0/0.917)*(sea_level+sl_change(t[i]))-Hght+10+(1-1.0/0.917)*Gamma_s[:,i][-1],xmin=0.8,linestyle='--',color='purple',linewidth=1.5,label=r'$h_f$')


        if j == 1:
            plt.plot(X/1000,Gamma_h[:,0]-Hght+ 10,color='k',linestyle='--',linewidth=1.5)
            plt.plot(X[Gamma_s[:,0]-bed(X)>tol]/1000,Gamma_s[:,0][Gamma_s[:,0]-bed(X)>tol],color='k',linestyle='--',linewidth=1.5)


        if j == 1 or j==3:
            plt.yticks([-2.5,0,2.5,10,15,20],[r'$-2.5$',r'$0$',r'$2.5$',r'$h_0$',r'$h_0+5$',r'$h_0+10$'],fontsize=16)
            plt.ylabel(r'$z$ (m)',fontsize=20)
        else:
           plt.yticks([-2.5,0,2.5,10,15,20],[r'$-2.5$',r'$0$',r'$2.5$',r'$h_0$',r'$h_0+5$',r'$h_0+10$'],fontsize=16)
           plt.gca().get_yaxis().set_ticklabels([])

        if j >2:
            plt.xlabel(r'$x$ (km)',fontsize=20)
            plt.xticks(2*np.array([0,1,2,3,4,5,6,7,8,9,10]),['-10','-8','-6','-4','-2','0','2','4','6','8','10'],fontsize=16)
        else:
            plt.xticks(2*np.array([0,1,2,3,4,5,6,7,8,9,10]),['','','','','','','','','','',''],fontsize=16)

        plt.yticks(fontsize=16)
        plt.ylim(np.min(bed(X))-2,20.0,8)
        plt.xlim(0,Lngth/1000)

        j+=1

    plt.tight_layout()
    lgd = plt.legend(fontsize=20,bbox_to_anchor=(0.95, -0.175),ncol=6)
    plt.savefig('fig3', bbox_extra_artists=(lgd,),bbox_inches='tight')
    plt.close()

#-------------------------------------------------------------------------------

def paperplot_fig4(s_mean,h_mean,lake_vol,x_right,x_left,dPw):

    t = np.linspace(0,t_final,num=int((t_final/3.154e7)*nt_per_year))

    B = 2*np.array([105,350,885,1020]) #Indices corresponding to free surface plots

    plt.figure(figsize=(13,6))
    plt.subplot(221)
    plt.annotate(r'(a)',xy=(-0.075,1.45),fontsize=20)
    plt.plot(t/3.154e7,Vol(t,1),'k',linewidth=3,label=r'true')
    plt.plot(t/3.154e7,lake_vol/lake_vol[0],'seagreen',linewidth=3,linestyle='--',label=r'computed')
    plt.plot(t[B]/3.154e7,Vol(t[B],1),'kD',markersize=10)
    plt.ylabel(r'$\mathcal{V}\,/\,\mathcal{V}_0$',fontsize=20)
    plt.gca().get_xaxis().set_ticklabels([])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(fontsize=16,loc='lower right',framealpha=1.0)


    plt.subplot(222)
    plt.annotate(r'(b)',xy=(-0.075,1.55),fontsize=20)
    p1, = plt.plot(t/3.154e7,s_mean-s_mean[0],'crimson',label=r'$\Delta\bar{s}$',linewidth=3)
    p2, = plt.plot(t/3.154e7,h_mean-h_mean[0],'royalblue',label=r'$\Delta\bar{h}$',linewidth=3)
    plt.ylabel(r'$z$ (m)',fontsize=20)
    plt.gca().get_xaxis().set_ticklabels([])
    plt.gca().yaxis.tick_right()
    plt.gca().yaxis.set_label_position("right")
    plt.ylim(-2,2)
    plt.xticks(fontsize=16)
    plt.yticks(np.arange(-2,2.1,1),fontsize=16)
    plt.legend(fontsize=16,loc='lower right',framealpha=1.0)

    plt.subplot(223)
    plt.annotate(r'(c)',xy=(-0.075,3.1),fontsize=20)
    plt.fill_between(t/(3.154e7),y1=x_left/1000.0 -5, y2=x_right/1000.0 - 5,facecolor='mediumslateblue',alpha=0.5,label='lake water')
    plt.plot(t/3.154e7,x_right/1000.0 - 5,'o',markersize=3,color='k')
    plt.plot(t/3.154e7,x_left/1000.0 -5 ,'o',markersize=3,color='k')
    plt.gca().ticklabel_format(axis='x',style='sci')
    plt.gca().yaxis.major.formatter._useMathText = True
    plt.xlabel(r'$t$ (yr)',fontsize=20)
    plt.ylabel(r'$x_\pm$ (km)',fontsize=20)
    plt.xticks(fontsize=16)
    plt.yticks(np.arange(-3.0,3.5,1.0),fontsize=16)
    plt.annotate(r'lake water',xy=(1.25,-0.25),fontsize=20)
    plt.ylim(-4,4)

    plt.subplot(224)
    plt.annotate(r'(d)',xy=(-0.075,0.075),fontsize=20)
    plt.axhline(y=0.0,linewidth=2,linestyle='-',color='navy')
    plt.plot(t/3.154e7,dPw,'ko',markersize=3)
    plt.gca().ticklabel_format(axis='x',style='sci')
    plt.gca().yaxis.major.formatter._useMathText = True
    plt.xlabel(r'$t$ (yr)',fontsize=20)
    plt.gca().yaxis.tick_right()
    plt.gca().yaxis.set_label_position("right")

    plt.ylabel(r'$\bar{p}_\mathrm{w}-p_\mathrm{o}$ (kPa)',fontsize=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.tight_layout()
    plt.savefig('fig4', bbox_inches='tight')
    plt.close()

#-------------------------------------------------------------------------------

def paperplot_fig5(Gamma_s,Gamma_h,x_left,x_right):
    X = X_fine

    # Sample time indices
    B = 2*np.array([105,350,885,1020])

    plt.figure(figsize=(12,8))
    j=1

    labels = [r'(a) lowstand',r'(b) filling',r'(c) highstand',r'(d) draining']

    for i in B:
        str1 = '22'+str(j)
        plt.subplot(int(str1))


        plt.annotate(labels[j-1],xy=(0.1,11.5),fontsize=20)
        plt.plot(X/1000,Gamma_h[:,i]-Hght+ 10,color='royalblue',linewidth=1.5,label=r'$h$')

        plt.fill_between(X/1000,y1=Gamma_s[:,i], y2=Gamma_h[:,i]-Hght+ 10,facecolor='aliceblue',alpha=1.0)
        plt.fill_between(X/1000,bed(X),Gamma_s[:,i],facecolor='slateblue',alpha=0.5)
        plt.fill_between(X/1000,-13*np.ones(np.size(X)),bed(X),facecolor='burlywood',alpha=1.0)

        plt.plot(X/1000,bed(X),color='k',linewidth=1.5,label=r'$\beta$')
        plt.plot(X[Gamma_s[:,i]-bed(X)>tol]/1000,Gamma_s[:,i][Gamma_s[:,i]-bed(X)>tol],color='crimson',linewidth=1.5,label=r'$s>\beta$')

        plt.plot(np.array([x_left[i]/1000]),np.array(np.min(bed(X))+1-2),marker='^',color='k',linestyle='None',markersize=10,label=r'$x_\pm$')
        plt.plot(np.array([x_right[i]/1000]),np.array(np.min(bed(X))+1-2),marker='^',markersize=10,color='k')

        if j == 1:
            plt.plot(X/1000,Gamma_h[:,0]-Hght+ 10,color='k',linestyle='--',linewidth=1.5)
            plt.plot(X[Gamma_s[:,0]-bed(X)>tol]/1000,Gamma_s[:,0][Gamma_s[:,0]-bed(X)>tol],color='k',linestyle='--',linewidth=1.5)



        if j == 1 or j==3:
            plt.yticks([-10,-5,0,7.5,10,12.5],[r'$-10$',r'$-5$',r'$0$',r'$h_0-2.5$',r'$h_0$',r'$h_0+2.5$'],fontsize=16)
            plt.ylabel(r'$z$ (m)',fontsize=20)
        else:
            plt.yticks([-10,-5,0,7.5,10,12.5],[r'',r'',r'',r'',r'',r''],fontsize=16)
            plt.gca().get_yaxis().set_ticklabels([])

        if j >2:
            plt.xlabel(r'$x$ (km)',fontsize=20)
            plt.xticks([0,1,2,3,4,5,6,7,8,9,10],['-5','-4','-3','-2','-1','0','1','2','3','4','5'],fontsize=16)
        else:
            plt.xticks([0,1,2,3,4,5,6,7,8,9,10],['','','','','','','','','','',''],fontsize=16)

        plt.yticks(fontsize=16)
        plt.ylim(np.min(bed(X))-2,14.0,8)
        plt.xlim(0,Lngth/1000)

        j+=1

    plt.tight_layout()
    lgd = plt.legend(fontsize=20,bbox_to_anchor=(0.55, -0.225),ncol=4)
    plt.savefig('fig5', bbox_extra_artists=(lgd,),bbox_inches='tight')
    #plt.show()
    plt.close()


#-------------------------------------------------------------------------------

def gabstract(Gamma_s,Gamma_h,x_left,x_right):
    X = X_fine

    # Sample time indices
    B = 2*np.array([105,350,1020,885])

    plt.figure(figsize=(12,10))
    j=1

    for i in B:
        str1 = '22'+str(j)
        plt.subplot(int(str1))

        plt.plot(X/1000,Gamma_h[:,i]-Hght+ 10,color='royalblue',linewidth=1.5,label=r'$h$')

        plt.fill_between(X/1000,y1=Gamma_s[:,i], y2=Gamma_h[:,i]-Hght+ 10,facecolor='aliceblue',alpha=1.0)
        plt.fill_between(X/1000,bed(X),Gamma_s[:,i],facecolor='slateblue',alpha=0.5)
        plt.fill_between(X/1000,-13*np.ones(np.size(X)),bed(X),facecolor='burlywood',alpha=1.0)

        plt.plot(X/1000,bed(X),color='k',linewidth=1.5,label=r'$\beta$')
        plt.plot(X[Gamma_s[:,i]-bed(X)>tol]/1000,Gamma_s[:,i][Gamma_s[:,i]-bed(X)>tol],color='crimson',linewidth=1.5,label=r'$s>\beta$')

        plt.plot(np.array([x_left[i]/1000]),np.array(np.min(bed(X))+1-2),marker='^',color='k',linestyle='None',markersize=10,label=r'$x_\pm$')
        plt.plot(np.array([x_right[i]/1000]),np.array(np.min(bed(X))+1-2),marker='^',markersize=10,color='k')

        plt.xticks([])
        plt.yticks([])

        plt.ylim(np.min(bed(X))-2,14.0,8)
        plt.xlim(0,Lngth/1000)

        j+=1

    plt.tight_layout()
    plt.savefig('gabstract')
    #plt.show()
    plt.close()

#-------------------------------------------------------------------------------
