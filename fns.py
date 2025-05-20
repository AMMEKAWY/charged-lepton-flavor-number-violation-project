import ROOT 
import warnings
from ROOT import TLorentzVector, TTree, TCanvas, TFile, TEfficiency
warnings.filterwarnings("ignore")

e_mass=0.00051 #GeV/c^2
mu_mass=0.10566 #GeV/c^2

f1=TFile("/afs/cern.ch/user/a/amekawy/summ/W_tau_eee/output_nanoaod_1.root")
f2=TFile("/afs/cern.ch/user/a/amekawy/summ/W_tau_mee/output_nanoaod_1.root")
f3=TFile("/afs/cern.ch/user/a/amekawy/summ/W_tau_mme/output_nanoaod_1.root")

pathes =[  "/afs/cern.ch/user/a/amekawy/summ/W_tau_eee/",
	 "/afs/cern.ch/user/a/amekawy/summ/W_tau_mee/", 
	 "/afs/cern.ch/user/a/amekawy/summ/W_tau_mme/"]
	 
channel=["t -> eee", "t -> mee", "t -> mme"]

tree1=f1.Get("Events;1")
tree2=f2.Get("Events;1")
tree3=f3.Get("Events;1")

branch_name = ["nGenPart","GenPart_pt","GenPart_eta","GenPart_phi"]
mu_rec_branch_name=["nMuon", "Muon_pt", "Muon_eta", "Muon_phi"]
e_rec_branch_name=["nElectron", "Electron_pt", "Electron_eta", "Electron_phi"]

N1=tree1.GetEntries()
N2=tree2.GetEntries()
N3=tree3.GetEntries()

trees=[tree1, tree2, tree3]
Ns=[N1, N2, N3]

#particleId1 e = 11, index=0
#particleId2 mu = 13, index=1
#parentId t = 15, index=2
#NoE == Number of Entries
#pin == particle index 

particleId=[11, 13, 15]

#=====================================================================

def plotter(firsthist, secondhist, nameA, nameB, plotname, path):
	
	global c1
	c1=ROOT.TCanvas("c1", plotname)
	c1.cd()
	firsthist.SetLineColor(ROOT.kRed)
	secondhist.SetLineColor(ROOT.kBlack)
	firsthist.SetLineStyle(1)  # solid line
	secondhist.SetLineStyle(2)  # dashed line
	firsthist.Draw()
	secondhist.Draw("same")

	legend = ROOT.TLegend(0.7, 0.8, 0.9, 0.9)
	legend.AddEntry(firsthist, nameA, "l")
	legend.AddEntry(secondhist, nameB, "l")
	legend.Draw()

	c1.SaveAs(path+plotname)
	pass

#=====================================================================

def effi(hist1, hist2, lable, path, filename="Efficiency.png"):

	c2=ROOT.TCanvas("c2", "Efficiency vs lepton p_{T}", 200, 10, 700, 500)
	eff=TEfficiency(hist1,hist2)
	eff.SetTitle(lable)
	eff.Draw()
	c2.SaveAs(path+filename)
	
#=====================================================================
#to avoid resetting the histograms in every iteration
def histcreator(name, name2, binnum, minimum, maximum):

	hist=ROOT.TH1F(name, name2, binnum, minimum, maximum)
	return hist

#=====================================================================

def finder(ttr, bbr_name, bbr_namePt, bbr_nameEta, bbr_namePhi, parentId, particleId1, particleId2 , NoE, genidx, hist_mu_pt, hist_mu_eta, hist_mu_phi, hist_e_pt, hist_e_eta, hist_e_phi , num=3):
	
	lv=TLorentzVector()
	ak=[]
	
	if abs(ttr.GenPart_pdgId[genidx]) in particleId and ttr.GenPart_status[genidx] == 1:
		if abs(ttr.GenPart_pdgId[ttr.GenPart_genPartIdxMother[genidx]]) == parentId:
			

			if abs(ttr.GenPart_pdgId[genidx])==particleId[0]:
			
			
				#electron		
				ak.append(11)
				b=e_rec_branch_name[0]
				lv.SetPtEtaPhiM(getattr(ttr, bbr_namePt)[genidx], getattr(ttr, bbr_nameEta)[genidx], getattr(ttr, bbr_namePhi)[genidx], e_mass)

				if abs(lv.Eta()) < 2.5:
					hist_e_pt.Fill(getattr(ttr, bbr_namePt)[genidx])	
					hist_e_eta.Fill(getattr(ttr, bbr_nameEta)[genidx])
					hist_e_phi.Fill(getattr(ttr, bbr_namePhi)[genidx])
					#h_accept_den_e_pt.Fill(getattr(ttr, bbr_namePt)[genidx])
					#h_accept_den_pt.Fill(getattr(ttr, bbr_namePt)[genidx])
					#print(ak)
			elif abs(ttr.GenPart_pdgId[genidx])==particleId[1]:
								
				#muon		
				ak.append(13)
				#print(ak)
				b=mu_rec_branch_name[0]
				lv.SetPtEtaPhiM(getattr(ttr, bbr_namePt)[genidx], getattr(ttr, bbr_nameEta)[genidx], getattr(ttr, bbr_namePhi)[genidx], mu_mass)
				
				if abs(lv.Eta()) < 2.5:
					#dc_mu_pt.append(lv)
					hist_mu_pt.Fill(getattr(ttr, bbr_namePt)[genidx])	
					hist_mu_eta.Fill(getattr(ttr, bbr_nameEta)[genidx])
					hist_mu_phi.Fill(getattr(ttr, bbr_namePhi)[genidx])	
					#h_accept_den_pt.Fill(getattr(ttr, bbr_namePt)[genidx])
					#print(lv.M())
	
	return lv, ak

#===============================================================

def matcher(partlist, tree, branch, branchpt, brancheta, branchphi, m):
	
	recomu = TLorentzVector()	
	rs=[]
	min_dR = 9999
	best_recidx = -1
	#for klm in range(particlex):
	for i in range(getattr(tree, branch)):
		
		recomu.SetPtEtaPhiM(getattr(tree, branchpt)[i],getattr(tree, brancheta)[i],getattr(tree, branchphi)[i],m)
		dR = partlist.DrEtaPhi(recomu)
		
		if dR<min_dR:
			min_dR = dR
			best_recidx = i
	
	return best_recidx, min_dR, recomu	

#=============================================================		

def arrs(N, ttr, bn=branch_name[0], num=3):	

	global mu_pt, mu_eta, mu_phi, e_pt, e_eta, e_phi

	mu_pt=ROOT.TH1F("Mu_pt", "muon_pt", 100, 0, 100)
	mu_eta=ROOT.TH1F("Mu_eta", "muon_eta",100, -10, 10)
	mu_phi=ROOT.TH1F("Mu_phi", "muon_phi", 100, -10, 10)

	e_pt=ROOT.TH1F("e_pt", "e_pt", 100, 0, 100)	
	e_eta=ROOT.TH1F("e_eta", "e_eta",100, -10, 10)
	e_phi=ROOT.TH1F("e_phi", "e_phi", 100, -10, 10)
		
	rec_mu_pt=ROOT.TH1F("rec_mu_pt", "muon_pt", 100, 0, 100)
	rec_mu_eta=ROOT.TH1F("rec_mu_eta", "muon_eta",100, -10, 10)
	rec_mu_phi=ROOT.TH1F("rec_mu_phi", "muon_phi", 100, -10, 10)

	rec_e_pt=ROOT.TH1F("rec_e_pt", "e_pt", 100, 0, 100)
	rec_e_eta=ROOT.TH1F("rec_e_eta", "e_eta", 100, -10, 10)
	rec_e_phi=ROOT.TH1F("rec_e_phi", "e_phi", 100, -10, 10)

	h_accept_num_pt = ROOT.TH1F('h_accept_num_pt', "acceptance efficiency numerator", 100, 0, 100)	
	h_accept_den_pt = ROOT.TH1F('h_accept_den_pt', "acceptance efficiency denominator", 100, 0, 100)	
		
	for i in range (N):

		ttr.GetEntry(i)
		dec=[]
		rec=[]
		
		for genid in range(getattr(ttr, bn)):
	
			ch, ip=finder(ttr, branch_name[0], branch_name[1],branch_name[2],branch_name[3], particleId[2], particleId[0], particleId[1], N, genid, mu_pt, mu_eta, mu_phi, e_pt, e_eta, e_phi)			
			
			if ch.Mag() == 0:
				continue
	
			if  len(ip) == 0:
				continue

			#print(abs(ch.Eta()))
			h_accept_den_pt.Fill(ch.Pt())
			
			if ip[0] == particleId[0]:
			
				index, r, reconst = matcher(ch, ttr, e_rec_branch_name[0], e_rec_branch_name[1], e_rec_branch_name[2], e_rec_branch_name[3], e_mass)	
				if r < 0.01:			
					rec_e_pt.Fill(getattr(ttr, e_rec_branch_name[1])[index])
					#print(getattr(ttr, e_rec_branch_name[2])[index])
					rec_e_eta.Fill(getattr(ttr, e_rec_branch_name[2])[index])
					rec_e_phi.Fill(getattr(ttr, e_rec_branch_name[3])[index])
					rec.append(reconst)
					#print(index, r)						
					#h_accept_num_e_pt.Fill(ch.Pt())
					h_accept_num_pt.Fill(ch.Pt())
			
			elif ip[0] == particleId[1]:	

				index, r, reconst = matcher(ch, ttr, mu_rec_branch_name[0], mu_rec_branch_name[1], mu_rec_branch_name[2], mu_rec_branch_name[3], mu_mass)	
				
				if r <0.01:
					rec_mu_pt.Fill(getattr(ttr, mu_rec_branch_name[1])[index])
					rec_mu_eta.Fill(getattr(ttr, mu_rec_branch_name[2])[index])
					rec_mu_phi.Fill(getattr(ttr, mu_rec_branch_name[3])[index])
					rec.append(reconst)				
					#print(index, r)
					#h_accept_num_mu_pt.Fill(ch.Pt())	
					h_accept_num_pt.Fill(ch.Pt())

	return rec_e_pt, rec_e_eta, rec_e_phi, rec_mu_pt, rec_mu_eta, rec_mu_phi, h_accept_num_pt, h_accept_den_pt 
	
#=============================================================	

def main():

	for i in range(len(trees)):

		print("=================", channel[i], "=================")

		rept, reeta, rephi, remupt, remueta, remuphi, numpt, denpt = arrs(Ns[i], trees[i])

		#plotter(firsthist (gen), secondhist (rec), nameA, nameB, plotname, path)
	
		if i == 0:
	
			plotter(e_eta, reeta, "e_eta", "rec_e_eta" , "e_eta.png", pathes[i])
			plotter(e_pt, rept, "e_pt", "rec_e_pt", "e_pt.png", pathes[i])
			plotter(e_phi, rephi, "e_phi", "rec_e_phi" , "e_phi.png", pathes[i])

		else:

			plotter(e_eta, reeta, "e_eta", "rec_e_eta" , "e_eta.png", pathes[i])
			plotter(e_pt, rept, "e_pt", "rec_e_pt", "e_pt.png", pathes[i])
			plotter(e_phi, rephi, "e_phi", "rec_e_phi" , "e_phi.png", pathes[i])
			
			plotter(mu_pt, remupt, "mu_pt", "rec_mu_pt", "mu_pt.png", pathes[i])
			plotter(mu_phi, remuphi, "mu_phi", "rec_mu_phi", "mu_phi.png", pathes[i])
			plotter(mu_eta, remueta, "mu_eta", "rec_mu_eta", "mu_eta.png", pathes[i])


		effi(numpt, denpt, "eff", pathes[i], filename="Efficiency.png")

	print("=================", "done", "=================")

#=============================================================	

if __name__ == "__main__":
	main()
	
	
