package team.ptm.model.service;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.LinkedHashMap;

import org.apache.http.client.HttpClient;
import org.apache.http.client.ResponseHandler;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.BasicResponseHandler;
import org.apache.http.impl.client.DefaultHttpClient;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.select.Elements;



public class FunctionMapping {

	public static void main(String[] args) {
		String root ="E:\\WM\\GOandKEGG\\data\\GO\\term\\";
		String spename = "HUMAN";
		FunctionalFeaturesUtil_1 getnm=new FunctionalFeaturesUtil_1();
		getnm.test(spename,getnm,0.05);
		System.out.println(root);
		try {
			GOMapping(root + "GoC.txt", 
					root + "mapping\\mapping_GoC.txt");
			GOMapping(root + "GoP.txt", 
					root + "mapping\\mapping_GoP.txt");
			GOMapping(root + "GoF.txt", 
					root + "mapping\\mapping_GoF.txt");
			KEGG_Mapping(getnm,spename,root + "KEGG.txt",
					root + "mapping\\mapping_KEGG.txt");
		} catch (Exception e) {
				e.printStackTrace();
		}
	}
	
	public static void KEGG_Mapping(FunctionalFeaturesUtil_1 getnm,String spename, String infile, String outfile) throws Exception {
		
		System.out.println(outfile);
		HttpClient httpclient = new DefaultHttpClient();
		BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
		bw.write("ID" + "\t" + "count"+"\t" +"M"+"\t"+"Pvalue"+ "\t"+"Term"+ "\t"+"n" + "\t"+"m"+ "\n");
		String tmpstr = "";
		while((tmpstr = br.readLine()) != null) {
			tmpstr = tmpstr.trim();
			if(tmpstr.equals(""))
				continue;
			String[] arr = tmpstr.split("\\s+");
			HttpGet submitpost = new HttpGet("http://www.genome.jp/dbget-bin/www_bget?pathway+" + arr[0].substring(5));
			ResponseHandler<String> handler = new BasicResponseHandler();
			String resstr = httpclient.execute(submitpost, handler);
			Document doc = Jsoup.parse(resstr.substring(resstr.indexOf("<nobr>Name</nobr></th>") + 22));
//			System.out.println(doc);
			Elements links = doc.select("div");	
//			System.out.println(links.get(0).text());
			bw.write(tmpstr + "\t" + links.get(0).text() + "\t"+getnm.getn + "\t"+getnm.getm+ "\n");
			bw.flush();
			System.out.println(arr[0] + "\t" + links.get(0).text());
	        submitpost.releaseConnection();
		}
		bw.flush();
		bw.close();
		br.close();
	}
	public static void KEGGmapping(String infile, String mappingfile, String outfile) throws Exception {
		LinkedHashMap<String, String> mappingMap = new LinkedHashMap<String, String>();
		BufferedReader br = new BufferedReader(new FileReader(new File(mappingfile)));
		String tmpstr = "";
		while((tmpstr = br.readLine()) != null) {
			tmpstr = tmpstr.trim();
			if(tmpstr.equals(""))
				continue;
			String[] arr = tmpstr.split("\t");
			System.out.println(arr[0] + "\t" + arr[1]);
			mappingMap.put(arr[0], arr[1].substring(0, arr[1].lastIndexOf('-') - 1));
		}
				
		br = new BufferedReader(new FileReader(new File(infile)));
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
		tmpstr = "";
		while((tmpstr = br.readLine()) != null) {
			tmpstr = tmpstr.trim();
			if(tmpstr.equals(""))
				continue;
			String[] arr = tmpstr.split("\\s+");
			bw.write(tmpstr + "\t" + mappingMap.get(arr[0]) + "\n");
		}
		bw.flush();
		bw.close();
	}
	
	
	public static void GOMapping(String infile, String outfile) throws Exception {
		System.out.println(outfile);
		HttpClient httpclient = new DefaultHttpClient();
		//HttpClient httpclient = new HttpClient();
		BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
		String tmpstr = "";
		while((tmpstr = br.readLine()) != null) {
			tmpstr = tmpstr.trim();
			if(tmpstr.equals(""))
				continue;
			String[] arr = tmpstr.split("\\s+");
			System.out.println(arr[0]);
			HttpGet submitpost = new HttpGet("http://amigo.geneontology.org/amigo/term/" + arr[0]);
			ResponseHandler<String> handler = new BasicResponseHandler();
			String resstr = httpclient.execute(submitpost, handler);
			Document doc = Jsoup.parse(resstr);
//			System.out.println(doc);
			Elements links = doc.select("h1");	
//			System.out.println(links.get(0).text());
			bw.write(tmpstr + "\t" + links.get(0).text() + "\n");
			bw.flush();
			System.out.println(arr[0] + "\t" + links.get(0).text());
	        submitpost.releaseConnection();
		}
		bw.flush();
		bw.close();
		br.close();
	}
	
	public static void InterProMapping(String infile, String outfile) throws Exception {
		HttpClient httpclient = new DefaultHttpClient();
		BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
		String tmpstr = "";
		while((tmpstr = br.readLine()) != null) {
			tmpstr = tmpstr.trim();
			if(tmpstr.equals(""))
				continue;
			String[] arr = tmpstr.split("\\s+");
			
			HttpGet submitpost = new HttpGet("http://www.ebi.ac.uk/interpro/entry/" + arr[0]);
			ResponseHandler<String> handler = new BasicResponseHandler();
			String resstr = httpclient.execute(submitpost, handler);
			Document doc = Jsoup.parse(resstr);
			Elements links = doc.select("h2.strapline");
			String retstr = links.get(0).text();
			bw.write(tmpstr + "\t" + retstr.substring(0, retstr.lastIndexOf('(') - 1) + "\n");
			bw.flush();			
			System.out.println(arr[0] + "\t" + retstr.substring(0, retstr.lastIndexOf('(') - 1));
	        submitpost.releaseConnection();
		}
		bw.flush();
		bw.close();
		br.close();
	}
	
	public static void PfamMapping(String infile, String outfile) throws Exception {
		HttpClient httpclient = new DefaultHttpClient();
		BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
		String tmpstr = "";
		while((tmpstr = br.readLine()) != null) {
			tmpstr = tmpstr.trim();
			if(tmpstr.equals(""))
				continue;
			String[] arr = tmpstr.split("\\s+");
			HttpGet submitpost = new HttpGet("http://pfam.xfam.org/family/" + arr[0]);
			ResponseHandler<String> handler = new BasicResponseHandler();
			String resstr = httpclient.execute(submitpost, handler);
			Document doc = Jsoup.parse(resstr);
			Elements links = doc.select("div#familySummaryBlock > div:first-child > h1");
			String retstr = links.get(0).text();
			bw.write(tmpstr + "\t" + retstr.substring(retstr.indexOf(":") + 2) + "\n");
			bw.flush();			
			System.out.println(arr[0] + "\t" + retstr.substring(retstr.indexOf(":") + 2));
	        submitpost.releaseConnection();
		}
		bw.flush();
		bw.close();
		br.close();
	}
	public static void PPIMapping(String infile, String outfile) throws Exception {
		HttpClient httpclient = new DefaultHttpClient();
		BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
		String tmpstr = "";
		while((tmpstr = br.readLine()) != null) {
			tmpstr = tmpstr.trim();
			if(tmpstr.equals(""))
				continue;
			String[] arr = tmpstr.split("\\s+");
			String protname = arr[0].split("\\.")[1];
			HttpGet submitpost = new HttpGet("http://www.uniprot.org/uniprot/?sort=score&query=" + protname);
			ResponseHandler<String> handler = new BasicResponseHandler();
			String resstr = httpclient.execute(submitpost, handler);
			Document doc= Jsoup.parse(resstr);
			Elements links = doc.select("div.protein_names > div.long");
			String retstr = "";
			try {
				retstr = links.get(0).text();
			} catch (IndexOutOfBoundsException e) {
				e.printStackTrace();
			}
			bw.write(tmpstr + "\t" + retstr + "\n");
			bw.flush();
			System.out.println(arr[0] + "\t" + retstr);
	        submitpost.releaseConnection();
		}
		bw.flush();bw.close();
		br.close();
	}
}
