#Samoel Giehl, Elias Nasr Naim Elias, Luciene. S. Delazari, Regiane Dalazoana
#samoelgiehl@gmail.com, elias_naim2008@hotmail.com, luciene@ufpr.br, regiane@ufpr.br
#Programa de Pós-Graduação em Ciências Geodésicas
#Universidade Federal do Paraná

import os
import wx
import geopandas as gpd
import numpy as np

wi = "Shapefile (*.shp)|*.shp|" \
            "All files (*.*)|*.*"

class Interface(wx.Frame):
    def __init__(self,parent):
        wx.Frame.__init__(self, parent, wx.ID_ANY, "TM_Proj",size=(450,300))
        self.currentDirectory = os.getcwd()        

        #Inserir texto na interface----------------------------------------------------------------------------------------

        wx.StaticText(self, 1, 'Projeção', (100,80))
        wx.StaticText(self, 1, 'Fuso (UTM)', (290,80))
        wx.StaticText(self, 1, 'Universidade Federal do Paraná', (120,5))
        wx.StaticText(self, 1, 'Programa de Pós-graduação em Ciências Geodésicas', (40,25))
                
        #Inserir o botão "Abrir Shapefile e Computar" --------------------------------------------------------------------------

        compute = wx.Button(self, label='Abrir Shapefile e Computar',size=(200,45),pos=(200,200))
        compute.Bind(wx.EVT_BUTTON, self.onOpenFile)

        #Inserir o botão "Fechar"--------------------------------------------------------------------------------------------
        
        close = wx.Button(self, label="Fechar",size=(120,45),pos=(30,200))
        close.Bind(wx.EVT_BUTTON, self.OnClose)
        
        #Inserir uma caixa de listagem de seleção múltipla-------------------------------------------------------------------

        lista1 = ["AUTO", "UTM", "LTM","RTM"]
        self.pro = wx.ComboBox(self,-1, "Escolha uma Projeção",size=(200,30),pos=(20,100))
        self.pro.SetItems(lista1)
     

        lista2 = ["MC -75° --> Fuso 18", "MC -69° --> Fuso 19", "MC -63° --> Fuso 20","MC -57° --> Fuso 21",
                  "MC -51° --> Fuso 22", "MC -45° --> Fuso 23","MC -39° --> Fuso 24","MC -33° --> Fuso 25"]
        self.fus = wx.ComboBox(self,-1, "Escolha o Fuso",size=(150,30),pos=(260,100))
        self.fus.SetItems(lista2)
    
    #Fechar programa--------------------------------------------------------------------------------------------------------
    
    def OnClose(self, event):
        self.Destroy()
 
    #Abrir shapefile e computar --------------------------------------------------------------------------------------------
    
    def onOpenFile(self, event):
        global proj,lon_mc,pr,deform
        dlg = wx.FileDialog(
            self, message="Escolha um Shapefile",
            defaultDir=self.currentDirectory, 
            defaultFile="",
            wildcard=wi,
            style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            paths = dlg.GetPaths()  
            for path in paths:
                pr= self.pro.GetValue()
                fuso=int(self.fus.GetValue()[17:19]) #Fuso de localização do shapefile
                sh = gpd.GeoDataFrame.from_file(path)
                long2=sh.bounds['maxx'][0]
                long1=sh.bounds['minx'][0]
                dif=np.abs(long2-long1)
                lat_cen=np.deg2rad(sh.iloc[0].geometry.centroid.y)
                lon_cen=np.deg2rad(sh.iloc[0].geometry.centroid.x)
                def round_to_odd(mmc):
                    mmc = int(mmc)
                    return mmc - 1 if mmc % 2 == 0 else mmc
                
                #Distorção de escala em função das coordenadas geodésicas (m)-----------------------------------------------
                
                def deformacao(lat_ce,dl):
                                a=6378137           #Semi-eixo maior do elipsoide GRS80
                                b=6356752.3141      #Semi-eixo menor do elipsoide GRS80
                                ee=((a**2-b**2)**(1/2))/b
                                n=ee**2*(np.cos(lat_cen))**2
                                t=np.tan(lat_cen)
                                m1=(dl**2)*(((np.cos(lat_cen))**2)/2)*(1+n**2)
                                m2=(dl**4)*(((np.cos(lat_cen))**4)/24)*(5-4*t**2+14*n**2+13*n**4-28*t**2*n**2+4*n**6-48*t**2*n**4-24*t**2*n**6)
                                m3=(dl**6)*(((np.cos(lat_cen))**6)/720)*(61-148*t**2+16*t**4)
                                m=1+m1+m2+m3
                                return m
                
                #Cálculo da longitude do MC----------------------------------------------------------------------------------
               
                lon_mc_utm=-(183-6*fuso)
                lon_mc_ltm=round((sh.iloc[0].geometry.centroid.x)*2)/2
                lon_mc_rtm=round_to_odd(sh.iloc[0].geometry.centroid.x)
                
                # Condições: LTM, RTM, UTM ou AUTO----------------------------------------------------------------------------
                
                if pr=='AUTO':
                        dlon_ltm=np.deg2rad(lon_mc_ltm)-lon_cen
                        dlon_rtm=np.deg2rad(lon_mc_rtm)-lon_cen
                        dlon_utm=np.deg2rad(lon_mc_utm)-lon_cen
                        deformm=np.array([abs(deformacao(lat_cen,dlon_ltm)*0.999995-1),abs(deformacao(lat_cen,dlon_rtm)*0.999995-1),abs(deformacao(lat_cen,dlon_utm)*0.9996-1)])
                        if dif<1:
                            mini=np.argmin(deformm)
                            if mini==0:
                                 pr='LTM'
                            if mini==1:
                                 pr='RTM'
                            elif mini==2:    
                                 pr='UTM'                         
                        elif dif>1 and dif<=2: 
                            mini=np.argmin(deformm[1:])
                            if mini==0:
                                 pr='RTM'
                            elif mini==1:    
                                 pr='UTM' 
                        else:
                            pr='UTM' 
                if pr=='UTM':
                        lon_mc=-(183-6*fuso)
                        dlon=np.deg2rad(lon_mc)-lon_cen
                        deform=deformacao(lat_cen,dlon)*0.9996
                        utm=("+proj=utm +south +zone=%s +datum=WGS84 +ellps=GRS80 +units=m +no_defs"%(fuso))
                        proj=sh.to_crs(utm)
                elif pr=='LTM':
                        lon_mc=round((sh.iloc[0].geometry.centroid.x)*2)/2
                        dlon=np.deg2rad(lon_mc)-lon_cen
                        deform=deformacao(lat_cen,dlon)*0.999995
                        ltm=('+proj=tmerc +lat_0=0 +lon_0=%s +k_0=0.999995 +x_0=200000 +y_0=5000000 +datum=WGS84 +ellps=GRS80 +units=m +no_defs' %(lon_mc))
                        proj=sh.to_crs(ltm)
                elif pr=='RTM':
                        lon_mc=round_to_odd(sh.iloc[0].geometry.centroid.x)
                        dlon=np.deg2rad(lon_mc)-lon_cen
                        deform=deformacao(lat_cen,dlon)*0.999995
                        rtm=('+proj=tmerc +lat_0=0 +lon_0=%s +k_0=0.999995 +x_0=400000 +y_0=5000000 +datum=WGS84 +ellps=GRS80 +units=m +no_defs' %(lon_mc))
                        proj=sh.to_crs(rtm)
                self.Destroy()
                self.Salve = Interface2(wx.Frame)
                self.Salve.Show()

#Abrir uma nova interface para salvar o arquivo de texto-----------------------------------------------------------------
     
class  Interface2(wx.Frame):
    def __init__(self,parent):
        super(Interface2, self).__init__(None, size=(400,200))
        self.SetTitle('finalizado com sucesso')
        compute = wx.Button(self, label="\tSalvar como\n\n\n Meridiano Central --> %s \n Projetado para --> %s \n k --> %s"
                            %(lon_mc,pr,deform),size=(80,25),pos=(0,0))
        compute.Bind(wx.EVT_BUTTON, self.OnSaveAs)
        self.currentDirectory = os.getcwd()
    def OnSaveAs(self,Interface2):
        dlg = wx.FileDialog(self, message="TM_Proj: Salvar...",defaultDir=self.currentDirectory, 
            defaultFile="", wildcard=wi, style=wx.FD_SAVE| wx.FD_OVERWRITE_PROMPT)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            proj.to_file(path, encoding='utf-8')
        dlg.Destroy()
        self.Destroy()
        
if __name__ == "__main__":
    app = wx.App(False)
    frame = Interface(None)
    frame.Show()
    app.MainLoop()