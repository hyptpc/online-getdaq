void test()
{
  const int N = 12;
  double s1[N];
  for(int i=0; i<N; i++) s1[i] = 0.5 + 0.5*i;
  double s2[N] = {0.28, 0.5, 0.74, 1, 1.24, 
    		  1.46, 1.7, 1.94, 2.18, 2.36, 
		  2.6, 2.78};
  TGraph* graph = new TGraph(N, s1, s2);
  TCanvas* c1 = new TCanvas("c1", "source", 800,600);
  c1->cd();
  graph->SetMarkerStyle(20);
  graph->SetMarkerColor(kOrange);
  graph->Draw("APL");
}
