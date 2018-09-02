package team.ptm.model.service;

/**
 * 功能：计算功能特征
 * 时间：2017-5-2
 */
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.ListIterator;
import jsc.distributions.Hypergeometric;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.io.FastaReaderHelper;

abstract class mapExtractHelper implements Runnable {

	LinkedHashMap<String, LinkedList<String>> retMap;
	HashMap<String, LinkedList<String>> tmpMap;
	LinkedHashMap<String, Integer> backgroundProtMap;
	BackgroundProt[] backgroundProtArray;
	String termsFile;
	private int samplesize;
	private int populationsize;
	private double pvalueThreshold; // 阈值

	/**
	 * @return the pvalueThreshold
	 */
	public double getPvalueThreshold() {
		return pvalueThreshold;
	}

	/**
	 * @param pvalueThreshold
	 *            the pvalueThreshold to set
	 */
	public void setPvalueThreshold(double pvalueThreshold) {
		this.pvalueThreshold = pvalueThreshold;
	}

	/**
	 * @return the samplesize
	 */
	public int getSamplesize() {
		return samplesize;
	}

	/**
	 * @return the populationsize
	 */
	public int getPopulationsize() {
		return populationsize;
	}

	/**
	 * @param populationsize
	 *            the populationsize to set
	 */
	public void setPopulationsize(int populationsize) {
		this.populationsize = populationsize;
	}

	/**
	 * @param backgroundProtMap
	 *            the backgroundProtMap to set
	 */
	public void setBackgroundProtMap(LinkedHashMap<String, Integer> backgroundProtMap) {
		this.backgroundProtMap = backgroundProtMap;
	}

	/**
	 * @param backgroundProtArray
	 *            the backgroundProtArray to set
	 */
	public void setBackgroundProtArray(BackgroundProt[] backgroundProtArray) {
		this.backgroundProtArray = new BackgroundProt[backgroundProtArray.length];
		for (int i = 0; i < backgroundProtArray.length; i++)
			this.backgroundProtArray[i] = backgroundProtArray[i].clone();
	}

	/**
	 * @param samplesize
	 *            the samplesize to set
	 */
	public void setSamplesize(int samplesize) {
		this.samplesize = samplesize;
	}

	/**
	 * @return the retMap
	 */
	public LinkedHashMap<String, LinkedList<String>> getRetMap() {
		return retMap;
	}

	/**
	 * @return the backgroundProtMap
	 */
	public LinkedHashMap<String, Integer> getBackgroundProtMap() {
		return backgroundProtMap;
	}

	/**
	 * @return the backgroundProtArray
	 */
	public BackgroundProt[] getBackgroundProtArray() {
		return backgroundProtArray;
	}

	public mapExtractHelper() {
		retMap = new LinkedHashMap<String, LinkedList<String>>();
		tmpMap = new HashMap<String, LinkedList<String>>();
	}

	public mapExtractHelper(String termsfile) {
		this();
		this.termsFile = termsfile;
	}

	public synchronized void addTerm(String uniprotid, String term) {
		if (tmpMap.containsKey(uniprotid)) {
			tmpMap.get(uniprotid).add(term);
			// System.out.println("加入Map中："+arrStr[0]+":"+arrStr[1]);
		} else {
			LinkedList<String> list = new LinkedList<String>();
			list.add(term);
			tmpMap.put(uniprotid, list);
			// System.out.println("加入Map中："+arrStr[0]+":"+arrStr[1]);
		}
	}

	public void run() {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(termsFile)));
			String tmpstr = "";
			while ((tmpstr = br.readLine()) != null) {
				tmpstr = tmpstr.trim();
				if (tmpstr == "")
					continue;
				String arr[] = tmpstr.split("\\s+");
				addTerm(arr[0], arr[1]);
			}
			br.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		Iterator<String> iter = backgroundProtMap.keySet().iterator();
		while (iter.hasNext()) {
			String key = iter.next();
			retMap.put(key, tmpMap.get(key));
		}
	}
}

class interproMapextractHelper extends mapExtractHelper implements Runnable {
	public interproMapextractHelper(String termsfile) {
		super(termsfile);
	}
}

class pfamMapextractHelper extends mapExtractHelper implements Runnable {
	public pfamMapextractHelper(String termsfile) {
		super(termsfile);
	}
}

class GoFMapextractHelper extends mapExtractHelper implements Runnable {
	public GoFMapextractHelper(String termsfile) {
		super(termsfile);
	}
}

class GoCMapextractHelper extends mapExtractHelper implements Runnable {
	public GoCMapextractHelper(String termsfile) {
		super(termsfile);
	}
}

class GoPMapextractHelper extends mapExtractHelper implements Runnable {
	public GoPMapextractHelper(String termsfile) {
		super(termsfile);
	}
}

class KEGGMapextractHelper extends mapExtractHelper implements Runnable {

	LinkedHashMap<String, String> mappingMap;
	String pathwayfile;
	String mappingFile;

	public KEGGMapextractHelper(String pathwayfile, String mappingfile) {
		super();
		this.pathwayfile = pathwayfile;
		this.mappingFile = mappingfile;
		mappingMap = new LinkedHashMap<String, String>();
	}

	public synchronized void addMap(String UniprotID, String KEGGID) {
		mappingMap.put(UniprotID, KEGGID);
	}

	public HashMap<String, LinkedList<String>> getPathway() {
		HashMap<String, LinkedList<String>> map = new HashMap<String, LinkedList<String>>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(pathwayfile));
			String tmpstr = "";
			while ((tmpstr = br.readLine()) != null) {
				if (tmpstr == "")
					continue;
				String[] arrStr = tmpstr.trim().split("\\s+");
				if (map.containsKey(arrStr[0])) {
					map.get(arrStr[0]).add(arrStr[1]);
				} else {
					LinkedList<String> list = new LinkedList<String>();
					list.add(arrStr[1]);
					map.put(arrStr[0], list);
				}
			}
			br.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return map;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Runnable#run()
	 */
	@Override
	public void run() {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(mappingFile)));
			String tmpstr = "";
			while ((tmpstr = br.readLine()) != null) {
				tmpstr = tmpstr.trim();
				if (tmpstr == "")
					continue;
				String arr[] = tmpstr.split("\\s+");
				addMap(arr[0], arr[1]);
			}
			br.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

		tmpMap = getPathway();
		Iterator<String> iter = backgroundProtMap.keySet().iterator();
		while (iter.hasNext()) {
			String key = iter.next();
			retMap.put(key, tmpMap.get(mappingMap.get(key)));
		}
	}
}

class ppiMapextractHelper extends mapExtractHelper implements Runnable {
	int ppiThreshold;
	String ppiFile;
	String mappingFile;

	public ppiMapextractHelper(int ppiThreshold, String ppiFile, String mappingFile) {
		super();
		this.ppiThreshold = ppiThreshold;
		this.ppiFile = ppiFile;
		this.mappingFile = mappingFile;
	}

	/*
	 * public Thread testandAdd(String filename,String UniprotID,String[]
	 * testArr) { //final HashMap<String, LinkedList<String>> map=this.tmpMap;
	 * return new Thread(); }
	 */
	/**
	 * 从Mappingfiles中得到所有的Uniprot-String映射文件
	 * 
	 * @param infile
	 * @return Mapping Files
	 */
	public HashMap<String, String> getMappingFiles() {
		class tempclass {
			String name;
			float percent;

			public tempclass(String nm, float perc) {
				this.name = nm;
				this.percent = perc;
			}
		}
		HashMap<String, tempclass> map = new HashMap<String, tempclass>();
		HashMap<String, String> retmap = new HashMap<String, String>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(mappingFile));
			br.readLine();
			String tmpstr = "";
			while ((tmpstr = br.readLine()) != null) {
				tmpstr = tmpstr.trim();
				if (tmpstr == "")
					continue;
				String[] arrStr = tmpstr.split("\\s+");
				if (map.containsKey(arrStr[0])) {
					if (map.get(arrStr[0]).percent < Float.parseFloat(arrStr[2])) {
						map.put(arrStr[0], new tempclass(arrStr[1], Float.parseFloat(arrStr[2])));
					}
				} else {
					map.put(arrStr[0], new tempclass(arrStr[1], Float.parseFloat(arrStr[2])));
				}
			}
			// BufferedWriter bw = new BufferedWriter(new FileWriter(new
			// File("D:/works/testMap.txt")));
			for (Iterator<String> iter = map.keySet().iterator(); iter.hasNext();) {
				String key = iter.next();
				tempclass tmp = map.get(key);
				retmap.put(key, tmp.name);
				// System.out.println(key+"\t"+tmp.name);
				// bw.write(key+"\t"+tmp.name+"\n");
			}
			// bw.flush();
			// bw.close();
			br.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return retmap;
	}

	/**
	 * 获取蛋白蛋白相互作用
	 * 
	 * @return 蛋白蛋白相互作用映射
	 */
	public HashMap<String, LinkedList<String>> getPpiFile() {
		HashMap<String, LinkedList<String>> map = new HashMap<String, LinkedList<String>>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(ppiFile));
			String tmpstr = "";
			while ((tmpstr = br.readLine()) != null) {
				tmpstr = tmpstr.trim();
				if (tmpstr == "")
					continue;
				String[] arrStr = tmpstr.split("\\s+");
				if (Integer.parseInt(arrStr[2]) >= ppiThreshold) {
					if (map.containsKey(arrStr[0])) {
						map.get(arrStr[0]).add(arrStr[1]);
						// System.out.println("加入Map中："+arrStr[0]+":"+arrStr[1]);
					} else {
						LinkedList<String> list = new LinkedList<String>();
						list.add(arrStr[1]);
						map.put(arrStr[0], list);
						// System.out.println("加入Map中："+arrStr[0]+":"+arrStr[1]);
					}
				}
			}
			br.close();
			// BufferedWriter bw=new BufferedWriter(new FileWriter(new
			// File("D:/works/testMap.txt")));
			// LinkedList<String> tmplist=map.get("9606.ENSP00000420824");
			// System.out.println(tmplist.size());
			// bw.flush();
			// bw.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		// System.out.println(map.size());
		return map;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Runnable#run()
	 */
	@Override
	public void run() {
		HashMap<String, String> mappingmap = getMappingFiles();
		tmpMap = getPpiFile();
		Iterator<String> iter = backgroundProtMap.keySet().iterator();
		while (iter.hasNext()) {
			String key = iter.next();
			retMap.put(key, tmpMap.get(mappingmap.get(key)));
		}
	}
}

class BackgroundProt implements Cloneable {
	private int idx;
	private boolean isPos;
	private LinkedList<Integer> termlist;

	/**
	 * @return the idx
	 */
	public int getIdx() {
		return idx;
	}

	/**
	 * @param isPos
	 *            the isPos to set
	 */
	public void setPos(boolean isPos) {
		this.isPos = isPos;
	}

	/**
	 * @return the isPos
	 */
	public boolean getisPos() {
		return isPos;
	}

	public BackgroundProt(int index, boolean flag) {
		this.idx = index;
		this.termlist = new LinkedList<Integer>();
		this.isPos = flag;
	}

	public Integer[] gettermsArray() {
		return termlist.toArray(new Integer[0]);
	}

	public void insertTermId(Integer termid) {
		this.termlist.add(termid);
	}

	public synchronized BackgroundProt clone() {
		BackgroundProt prot = null;
		try {
			prot = (BackgroundProt) super.clone();
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		}
		if (prot.termlist != null) {
			LinkedList<Integer> tmplist = prot.termlist;
			Iterator<Integer> iter = tmplist.iterator();
			prot.termlist = new LinkedList<Integer>();
			while (iter.hasNext()) {
				prot.termlist.add(iter.next());
			}
		}
		return prot;
	}
}

class Term {
	private int idx;
	int countPos;
	int countBackground;
	private double pvalue;
	private double pPostive;
	private double pNegative;
	private LinkedList<Integer> UniprotAccID;

	/**
	 * @return the pvalue
	 */
	public double getPvalue() {
		return pvalue;
	}

	/**
	 * @return the pPostive
	 */
	public double getpPostive() {
		return pPostive;
	}

	/**
	 * @return the pNegative
	 */
	public double getpNegative() {
		return pNegative;
	}

	/**
	 * @param idx
	 *            the idx to set
	 */
	public void setIdx(int idx) {
		this.idx = idx;
	}

	/**
	 * @return the idx
	 */
	public int getIdx() {
		return idx;
	}

	/**
	 * @return the uniprotAccID
	 */
	public LinkedList<Integer> getUniprotAccID() {
		return UniprotAccID;
	}

	public int countUniprotAccList() {
		return UniprotAccID.size();
	}

	public Term() {
		UniprotAccID = new LinkedList<Integer>();
		countPos = 0;
		countBackground = 0;
	}

	public void insertUniprotAccID(Integer Accid) {
		UniprotAccID.add(Accid);
	}

	public double calculatePvalue(int samplesize, int populationsize, int hypersize) {
		Hypergeometric hyper = new Hypergeometric(samplesize, populationsize, countBackground);
		pvalue = hyper.pdf(countPos);
		return hypersize * pvalue;
		// return pvalue;
	}

	public void addcountPos() {
		countPos++;
	}

	public void addcountBackground() {
		countBackground++;
	}

	public void calculatep(int populationsize, int samplesize) {
		this.pPostive = Math.log((1.0 * countPos / samplesize) / (1.0 * countBackground / populationsize))
				/ Math.log(2.0);
		this.pNegative = Math.log((1 - (1.0 * countPos / samplesize)) / (1 - (1.0 * countBackground / populationsize)))
				/ Math.log(2.0);
	}
}

class GeneralFunctionalFeatures implements Runnable {
	HashMap<String, Integer> backgroundProtMap;
	BackgroundProt[] BackgroundProtArray;
	LinkedHashMap<String, LinkedList<String>> termslistMap;
	private int samplesize;
	private int populationsize;
	double pvalueThreshold;
	double[] Weight;
	double[] result;
	LinkedHashMap<String, Term> calculatedTermsMap;

	/**
	 * @return the result
	 */
	public double[] getResult() {
		return result;
	}

	public GeneralFunctionalFeatures(mapExtractHelper helper) {
		this.backgroundProtMap = helper.getBackgroundProtMap();
		BackgroundProtArray = helper.getBackgroundProtArray();
		this.termslistMap = helper.getRetMap();
		this.samplesize = helper.getSamplesize();
		this.populationsize = helper.getPopulationsize();
		this.pvalueThreshold = helper.getPvalueThreshold();
		result = new double[populationsize];
	}

	public void constructDataStructure() {
		LinkedHashMap<String, Term> termsMap = new LinkedHashMap<String, Term>();
		Iterator<String> iter = termslistMap.keySet().iterator();
		// 以下为带入超几何分布时只计算正样本蛋白，如需计算背景蛋白则删掉
		HashMap<String, String> hypercountMap = new HashMap<String, String>();
		// 结束
		while (iter.hasNext()) {
			String key = iter.next(); // 蛋白的key
			LinkedList<String> termsLinkedList;
			if ((termsLinkedList = termslistMap.get(key)) == null)
				continue;
			ListIterator<String> iter2 = termsLinkedList.listIterator();
			while (iter2.hasNext()) {
				String strterm = iter2.next();
				if (!termsMap.containsKey(strterm)) {
					termsMap.put(strterm, new Term());
				}
				Term tmpterm = termsMap.get(strterm);
				if (BackgroundProtArray[backgroundProtMap.get(key)].getisPos() == true) {
					tmpterm.addcountPos();
				}
				tmpterm.addcountBackground();
				tmpterm.insertUniprotAccID(BackgroundProtArray[backgroundProtMap.get(key)].getIdx());
				// 以下为带入超几何分布时只计算正样本蛋白，如需计算背景蛋白则替换下面一个注视行
				if (BackgroundProtArray[backgroundProtMap.get(key)].getisPos() == true)
					hypercountMap.put(strterm, " ");
				// 结束
			}
		}
		// 以下为带入超几何分布时只计算正样本蛋白，如需计算背景蛋白则替换下面一个注视行
		int hypersize = hypercountMap.size();
		// int hypersize=termsMap.size();
		// 结束
		System.out.println("Terms单独的Map(只是正样本)构造完毕,长度为：" + hypersize);
		Iterator<String> iter3 = termsMap.keySet().iterator();
		int count = 0;
		while (iter3.hasNext()) {
			String key = iter3.next();
			if ((termsMap.get(key).calculatePvalue(samplesize, populationsize, hypersize)) < pvalueThreshold) {
				// 测试System.out.print(key+"
				// "+termsMap.get(key).calculatePvalue(samplesize,
				// populationsize,hypersize)+" ");
				Term tmpTerm = termsMap.get(key);
				tmpTerm.setIdx(count);
				for (Iterator<Integer> iter2 = tmpTerm.getUniprotAccID().iterator(); iter2.hasNext();) {
					BackgroundProtArray[iter2.next()].insertTermId(count);
				}
				count++;
				tmpTerm.calculatep(populationsize, samplesize);
				// 测试System.out.println(samplesize+" "+populationsize+"
				// "+tmpTerm.countPos+" "+tmpTerm.countBackground+"
				// "+tmpTerm.getpNegative()+" "+tmpTerm.getpPostive());
			} else {
				iter3.remove();
			}
			// System.out.println(key+":"+p);
		}
		System.out.println("总共被选出来的term数量：" + termsMap.size());
		/*
		 * String key[]=termsMap.keySet().toArray(new String[0]); for(int
		 * i=0;i<key.length;i++) { System.out.println(key[i]); }
		 */
		calculatedTermsMap = termsMap;
	}

	public double[] calculateWeight() {
		int[][] Matrix;
		String key[] = calculatedTermsMap.keySet().toArray(new String[0]);
		Matrix = new int[this.calculatedTermsMap.size()][this.calculatedTermsMap.size()];
		for (int i = 0; i < key.length; i++) {
			for (int j = i + 1; j < key.length; j++)
				Matrix[i][j] = calculatedTermsMap.get(key[i]).countUniprotAccList()
						+ calculatedTermsMap.get(key[j]).countUniprotAccList();
		}
		for (int i = 0; i < key.length; i++) {
			for (int j = i + 1; j < key.length; j++)
				Matrix[j][i] = Matrix[i][j];
		}
		for (int k = 0; k < BackgroundProtArray.length; k++) {
			Integer[] array = BackgroundProtArray[k].gettermsArray();
			for (int i = 0; i < array.length; i++)
				for (int j = i + 1; j < array.length; j++) {
					Matrix[array[i]][array[j]]--;
					Matrix[array[j]][array[i]] -= 2;
				}
		}
		double[][] matrixDouble = new double[this.calculatedTermsMap.size()][this.calculatedTermsMap.size()];
		for (int i = 0; i < key.length; i++) {
			for (int j = i + 1; j < key.length; j++)
				matrixDouble[i][j] = 1.0 * Matrix[j][i] / Matrix[i][j];
		}
		double[] weight = new double[Matrix.length];
		for (int i = 0; i < matrixDouble.length; i++) {
			double sum = 0;
			for (int j = 0; j < i; j++) {
				sum += matrixDouble[j][i];
			}
			for (int j = i + 1; j < matrixDouble[i].length; j++) {
				sum += matrixDouble[i][j];
			}
			weight[i] = sum;
		}
		this.Weight = weight;
		// for(int i=0; i<weight.length; i++) {
		// //测试System.out.println(weight[i]);
		// }
		return this.Weight;
	}

	public double calculateoneScore(int idx) {
		Iterator<String> iter = calculatedTermsMap.keySet().iterator();
		String[] keys = calculatedTermsMap.keySet().toArray(new String[0]);
		double sum = 0;
		int count = 0;
		while (iter.hasNext()) {
			sum += 1.0 / (Weight[count++] + 1) * calculatedTermsMap.get(iter.next()).getpNegative();
		}
		Integer[] termsarr = BackgroundProtArray[idx].gettermsArray();
		for (int i = 0; i < termsarr.length; i++) {
			// 测试System.out.println("蛋白"+idx+"的"+termsarr[i]);
			Term term = calculatedTermsMap.get(keys[i]);
			// 测试System.out.println(term.getpPostive()+" "+term.getpNegative()+"
			// "+Weight[termsarr[i]]);
			sum = sum + 1.0 / (Weight[termsarr[i]] + 1) * (term.getpPostive() - term.getpNegative());
		}
		return sum;
	}

	public double calculateoneScore(String BackgroundProtkey) {
		return calculateoneScore(backgroundProtMap.get(BackgroundProtkey));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Runnable#run()
	 */
	@Override
	public void run() {
		constructDataStructure();
		calculateWeight();

		for (int i = 0; i < BackgroundProtArray.length; i++) {
			result[i] = calculateoneScore(i);
		}
		/*
		 * 以下为测试信息，获取矩阵，行为背景蛋白，列为term，输出文件testjaccard.txt int[][] arrayint=new
		 * int[86728][6]; for(int idx=0;idx<BackgroundProtArray.length;idx++) {
		 * Integer[] termsarr=BackgroundProtArray[idx].gettermsArray(); for(int
		 * i=0;i<termsarr.length;i++) {
		 * System.out.println("蛋白"+idx+"的"+termsarr[i]);
		 * arrayint[idx][termsarr[i]]=1; } } try { BufferedWriter bw = new
		 * BufferedWriter(new FileWriter(new File("testjaccard.txt"))); for(int
		 * i=0;i<arrayint.length;i++) { for(int j=0;j<arrayint[i].length;j++) {
		 * bw.write(arrayint[i][j]+" "); } bw.write("\n"); } bw.flush();
		 * bw.close(); } catch (Exception e) { e.printStackTrace(); }
		 */// 测试结束
	}
}

/**
 * 
 * @author
 *
 */
class FunctionalFeatures {
	String backgroundProteinFile;
	String postiveProtFile;
	double pvalueThreshold; // 阈值
	int samplesize;
	int populationsize;
	LinkedHashMap<String, Integer> backgroundProtMap;
	LinkedHashMap<String, Term> calculatedTermsList;
	LinkedHashMap<String, double[]> scores;
	LinkedHashMap<String, mapExtractHelper> jobs;
	BackgroundProt[] backgroundProtArray;

	
	public int mappingn;
	public int mappingm;
	/**
	 * @return the scores
	 */
	public LinkedHashMap<String, double[]> getScores() {
		return scores;
	}

	public LinkedHashMap<String, LinkedHashMap<String, Double>> getScoresbyHashMap() {
		Iterator<String> iter = scores.keySet().iterator();
		LinkedHashMap<String, LinkedHashMap<String, Double>> scoreMap = new LinkedHashMap<String, LinkedHashMap<String, Double>>();
		while (iter.hasNext()) {
			int count = 0;
			String key = iter.next();
			double[] score = scores.get(key);

			Iterator<String> iter2 = backgroundProtMap.keySet().iterator();
			LinkedHashMap<String, Double> protMap = new LinkedHashMap<String, Double>();
			while (iter2.hasNext()) {
				String backgroundkey = "";
				backgroundkey = iter2.next();
				// System.out.println(key+" "+backgroundkey+" "+score[count]);
				protMap.put(backgroundkey, score[count++]);
			}
			scoreMap.put(key, protMap);
		}
		return scoreMap;
	}

	public FunctionalFeatures(String backgroundProteinFilePath, String postiveProtFilePath, double pvalueThresh) {
		this.backgroundProteinFile = backgroundProteinFilePath;
		this.postiveProtFile = postiveProtFilePath;
		this.pvalueThreshold = pvalueThresh;
		this.backgroundProtMap = new LinkedHashMap<String, Integer>();
		jobs = new LinkedHashMap<String, mapExtractHelper>();
		scores = new LinkedHashMap<String, double[]>();
		prepare();
	}

	public void prepare() {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(backgroundProteinFile)));
			String tmpstr = "";
			int countBackground = 0;
			while ((tmpstr = br.readLine()) != null) {
				tmpstr = tmpstr.trim();
				if (tmpstr == "")
					continue;
				this.backgroundProtMap.put(tmpstr, countBackground++);
				// System.out.println("put"+tmpstr);
			}
			br.close();
			int countPos = 0;
			br = new BufferedReader(new FileReader(postiveProtFile));
			LinkedList<Integer> posProt = new LinkedList<Integer>();
			while ((tmpstr = br.readLine()) != null) {
				tmpstr = tmpstr.trim();
				if (tmpstr == "")
					continue;
				countPos++;
				if (!backgroundProtMap.containsKey(tmpstr)) {
					backgroundProtMap.put(tmpstr, countBackground++);
				}
				posProt.add(backgroundProtMap.get(tmpstr));
			}
			this.samplesize = countPos;
			this.populationsize = countBackground;
			backgroundProtArray = new BackgroundProt[countBackground];
			for (int i = 0; i < countBackground; i++) {
				backgroundProtArray[i] = new BackgroundProt(i, false);
			}
			for (Iterator<Integer> iter = posProt.iterator(); iter.hasNext();) {
				backgroundProtArray[iter.next()].setPos(true);
			}
			br.close();
		} catch (Exception e) {
			// PrintStackTrace: handle exception
			e.printStackTrace();
		}
		// 以下是测试信息，以后可以删掉
		mappingn=backgroundProtArray.length;
		System.out.println("总蛋白数量：" + mappingn);
		int countpos = 0;
		for (int i = 0; i < backgroundProtArray.length; i++) {
			if (backgroundProtArray[i].getisPos() == true) {
				countpos++;
			}
		}
		System.out.println("其中标注的正样本蛋白数量为：" + countpos);
		mappingm=countpos;
		// 测试信息结束
	}

	public void docalculate() {
		// 运行每个helper的run，完成整理工作
		System.out.println("开始运行每个helper的整理工作");
		Iterator<String> jobiter = jobs.keySet().iterator();
		Thread[] threads = new Thread[jobs.size()];
		int count = 0;
		while (jobiter.hasNext()) {
			threads[count] = new Thread(jobs.get(jobiter.next()));
			threads[count++].start();
		}
		for (int i = 0; i < threads.length; i++) {
			try {
				threads[i].join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		// 创建GeneralFunctionalFeatures进行计算
		System.out.println("开始进行计算");
		Iterator<String> jobiter2 = jobs.keySet().iterator();
		GeneralFunctionalFeatures[] features = new GeneralFunctionalFeatures[jobs.size()];
		count = 0;
		while (jobiter2.hasNext()) {
			mapExtractHelper helper = jobs.get(jobiter2.next());
			features[count] = new GeneralFunctionalFeatures(helper);
			threads[count] = new Thread(features[count]);
			threads[count].start();
			count++;
		}
		for (int i = 0; i < threads.length; i++) {
			try {
				threads[i].join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		// 获取结果
		System.out.println("开始获取结果");
		Iterator<String> jobiter3 = jobs.keySet().iterator();
		count = 0;
		while (jobiter3.hasNext()) {
			String key = jobiter3.next();
			// 输出超几何分布结果
			LinkedHashMap<String, Term> map = features[count].calculatedTermsMap;
			try {
				BufferedWriter bw = new BufferedWriter(
						new FileWriter(new File("E:\\WM\\GOandKEGG\\data\\GO\\term\\" + key + ".txt")));
				for (String termkey : map.keySet()) {
					Term term = map.get(termkey);
					bw.write(termkey + "\t" + term.countPos + "\t" + term.countBackground + "\t" + term.getPvalue()
							+ "\n");
				}
				bw.flush();
				bw.close();
			} catch (IOException e) {
				e.printStackTrace();
			}

			// 完成输出超几何分布结果
			scores.put(key, features[count++].result);
		}
	}

	public void addPfam(String Pfamtermsfile) {
		pfamMapextractHelper pfam = new pfamMapextractHelper(Pfamtermsfile);
		pfam.setBackgroundProtArray(backgroundProtArray);
		LinkedHashMap<String, Integer> tmpmap = new LinkedHashMap<String, Integer>();
		tmpmap.putAll(backgroundProtMap);
		pfam.setBackgroundProtMap(tmpmap);
		pfam.setSamplesize(samplesize);
		pfam.setPopulationsize(populationsize);
		pfam.setPvalueThreshold(this.pvalueThreshold);
		jobs.put("Pfam", pfam);
	}

	public void addInterpro(String Interprotermsfile) {
		interproMapextractHelper interpro = new interproMapextractHelper(Interprotermsfile);
		interpro.setBackgroundProtArray(backgroundProtArray);
		LinkedHashMap<String, Integer> tmpmap = new LinkedHashMap<String, Integer>();
		tmpmap.putAll(backgroundProtMap);
		interpro.setBackgroundProtMap(tmpmap);
		interpro.setSamplesize(samplesize);
		interpro.setPopulationsize(populationsize);
		interpro.setPvalueThreshold(this.pvalueThreshold);
		jobs.put("Interpro", interpro);
	}

	public void addGo(String GoCtermsfile, String GoFtermsfile, String GoPtermsfile) {
		GoCMapextractHelper goc = new GoCMapextractHelper(GoCtermsfile);
		goc.setBackgroundProtArray(backgroundProtArray);
		LinkedHashMap<String, Integer> tmpmap = new LinkedHashMap<String, Integer>();
		tmpmap.putAll(backgroundProtMap);
		goc.setBackgroundProtMap(tmpmap);
		goc.setSamplesize(samplesize);
		goc.setPopulationsize(populationsize);
		goc.setPvalueThreshold(this.pvalueThreshold);
		jobs.put("GoC", goc);
		GoFMapextractHelper gof = new GoFMapextractHelper(GoFtermsfile);
		gof.setBackgroundProtArray(backgroundProtArray);
		tmpmap = new LinkedHashMap<String, Integer>();
		tmpmap.putAll(backgroundProtMap);
		gof.setBackgroundProtMap(tmpmap);
		gof.setSamplesize(samplesize);
		gof.setPopulationsize(populationsize);
		gof.setPvalueThreshold(this.pvalueThreshold);
		jobs.put("GoF", gof);
		GoPMapextractHelper gop = new GoPMapextractHelper(GoPtermsfile);
		gop.setBackgroundProtArray(backgroundProtArray);
		tmpmap = new LinkedHashMap<String, Integer>();
		tmpmap.putAll(backgroundProtMap);
		gop.setBackgroundProtMap(tmpmap);
		gop.setSamplesize(samplesize);
		gop.setPopulationsize(populationsize);
		gop.setPvalueThreshold(this.pvalueThreshold);
		jobs.put("GoP", gop);
	}

	public void addKEGG(String pathwayfile, String KEGGMappingfile) {
		KEGGMapextractHelper kegg = new KEGGMapextractHelper(pathwayfile, KEGGMappingfile);
		kegg.setBackgroundProtArray(backgroundProtArray);
		LinkedHashMap<String, Integer> tmpmap = new LinkedHashMap<String, Integer>();
		tmpmap.putAll(backgroundProtMap);
		kegg.setBackgroundProtMap(tmpmap);
		kegg.setSamplesize(samplesize);
		kegg.setPopulationsize(populationsize);
		kegg.setPvalueThreshold(this.pvalueThreshold);
		jobs.put("KEGG", kegg);
	}

	public void addPpi(int ppiThreshold, String ppiFile, String mappingFile) {
		ppiMapextractHelper ppi = new ppiMapextractHelper(ppiThreshold, ppiFile, mappingFile);
		ppi.setBackgroundProtArray(backgroundProtArray);
		LinkedHashMap<String, Integer> tmpmap = new LinkedHashMap<String, Integer>();
		tmpmap.putAll(backgroundProtMap);
		ppi.setBackgroundProtMap(tmpmap);
		ppi.setSamplesize(samplesize);
		ppi.setPopulationsize(populationsize);
		System.out.println(samplesize + " " + populationsize);
		ppi.setPvalueThreshold(this.pvalueThreshold);
		jobs.put("PPI", ppi);
	}
}

class TXTExtractHelper {
	LinkedHashMap<String, Integer> backgroundProtMap;
	String txtFilePath;
	String outputPath;
	String Interproterms;
	String Pfamterms;
	String GoFterms;
	String GoCterms;
	String GoPterms;
	String KEGGMap;

	public TXTExtractHelper(LinkedHashMap<String, Integer> backgroundProtMap, String txtFilePath) {
		this.backgroundProtMap = backgroundProtMap;
		this.txtFilePath = txtFilePath;
		this.Interproterms = "";
		this.Pfamterms = "";
		this.GoFterms = "";
		this.GoCterms = "";
		this.GoPterms = "";
		this.KEGGMap = "";
	}

	public void extract() {
		// 读取txt文件的注释，并加入到相应的helper中
		System.out.println("开始读取txt文件");
		Iterator<String> keyiter = this.backgroundProtMap.keySet().iterator();
		while (keyiter.hasNext()) {
			String key = keyiter.next();
			try {
				BufferedReader br = new BufferedReader(
						new FileReader(new File(this.txtFilePath + File.separator + key + ".txt")));
				String tmpline = "";
				// System.out.println(key);
				while ((tmpline = br.readLine()) != null) {
					tmpline = tmpline.trim();
					String[] testArr = tmpline.split("\\s+");
					if (tmpline.startsWith("DR") && testArr.length >= 3) {
						if (testArr[1].equals("InterPro;")) {
							Interproterms += key + "\t" + testArr[2].substring(0, testArr[2].length() - 1) + "\n";
						} else if (testArr[1].equals("Pfam;")) {
							Pfamterms += key + "\t" + testArr[2].substring(0, testArr[2].length() - 1) + "\n";
						} else if (testArr[1].equals("GO;") && testArr[3].startsWith("F:")) {
							GoFterms += key + "\t" + testArr[2].substring(0, testArr[2].length() - 1) + "\n";
						} else if (testArr[1].equals("GO;") && testArr[3].startsWith("C:")) {
							GoCterms += key + "\t" + testArr[2].substring(0, testArr[2].length() - 1) + "\n";
						} else if (testArr[1].equals("GO;") && testArr[3].startsWith("P:")) {
							GoPterms += key + "\t" + testArr[2].substring(0, testArr[2].length() - 1) + "\n";
						} else if (testArr[1].equals("KEGG;")) {
							KEGGMap += key + "\t" + testArr[2].substring(0, testArr[2].length() - 1) + "\n";
						}
					}
				}
				br.close();
			} catch (Exception e) {
				System.out.println("找不到"+key+".txt");
				//e.printStackTrace();
			}
		}
	}

	public void extractandOut(String outputPath) {
		extract();
		try {
			BufferedWriter fw = new BufferedWriter(
					new FileWriter(new File(outputPath + File.separator + "interproTerms.txt")));
			fw.write(Interproterms);
			fw.flush();
			fw.close();
			fw = new BufferedWriter(new FileWriter(new File(outputPath + File.separator + "pfamTerms.txt")));
			fw.write(Pfamterms);
			fw.flush();
			fw.close();
			fw = new BufferedWriter(new FileWriter(new File(outputPath + File.separator + "GoCTerms.txt")));
			fw.write(GoCterms);
			fw.flush();
			fw.close();
			fw = new BufferedWriter(new FileWriter(new File(outputPath + File.separator + "GoFTerms.txt")));
			fw.write(GoFterms);
			fw.flush();
			fw.close();
			fw = new BufferedWriter(new FileWriter(new File(outputPath + File.separator + "GoPTerms.txt")));
			fw.write(GoPterms);
			fw.flush();
			fw.close();
			fw = new BufferedWriter(new FileWriter(new File(outputPath + File.separator + "KEGGMap.txt")));
			fw.write(KEGGMap);
			fw.flush();
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}

/**
 * @author Tom
 *
 */
public class FunctionalFeaturesUtil_1 {
	/**
	 * 从全部蛋白蛋白相互作用文件中获取某一物种的蛋白蛋白相互作用
	 * 
	 * @param LinkFile
	 *            全部蛋白蛋白相互作用文件
	 * @param outfile
	 *            输出文件
	 * @param specieNo
	 *            物种编号
	 */
	public int getm;
	public int getn;
	public static void getLinksOfOneSpecie(String LinkFile, String outfile, String specieNo) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(LinkFile)));
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
			String tmpstr = "";
			while ((tmpstr = br.readLine()) != null) {
				tmpstr = tmpstr.trim();
				if (tmpstr == "")
					continue;
				if (tmpstr.startsWith(specieNo)) {
					bw.write(tmpstr + "\n");
				}
			}
			bw.flush();
			bw.close();
			br.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static String getoneFileUniprot2String(String specie, String specieNo, String Acc, String txtFilepath) {
		String target = "";
		try {
			BufferedReader br = new BufferedReader(
					new FileReader(new File(txtFilepath + File.separator + Acc + ".txt")));
			String tmpstr = "";
			if (specie == "Salmonella" || specie == "Yeast") {
				while ((tmpstr = br.readLine()) != null) {
					tmpstr = tmpstr.trim();
					if (tmpstr == "")
						continue;
					if (tmpstr.startsWith("GN ")) {
						int Wordidx = tmpstr.indexOf("OrderedLocusNames");
						if (Wordidx == -1)
							continue;
						target = tmpstr.substring(tmpstr.indexOf("=", Wordidx) + 1, tmpstr.indexOf(";", Wordidx));
						break;
					}
				}

			}
			if (specie == "Ecoli") {
				while ((tmpstr = br.readLine()) != null) {
					tmpstr = tmpstr.trim();
					if (tmpstr == "")
						continue;
					if (tmpstr.startsWith("GN ")) {
						int Wordidx = tmpstr.indexOf("OrderedLocusNames");
						if (Wordidx == -1)
							continue;
						target = tmpstr.substring(tmpstr.indexOf("=", Wordidx) + 1, tmpstr.indexOf(";", Wordidx));
						String[] arr = target.split(", ");
						for (int i = 0; i < arr.length; i++) {
							target = "";
							if (arr[i].matches("b\\d+")) {
								target = arr[i];
								break;
							}
						}
						break;
					}
				}
			}
			if (specie == "Rat") {
				while ((tmpstr = br.readLine()) != null) {
					tmpstr = tmpstr.trim();
					if (tmpstr == "")
						continue;
					if (tmpstr.startsWith("DR   Ensembl")) {
						target = tmpstr.split("; ")[2];
						target = target.trim();
						if (target == "")
							System.out.println(tmpstr);
						break;
					}
				}
			}
			br.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.err.println(Acc);
		}
		if (target == "")
			return "";
		else
			return Acc + "\t" + specieNo + "." + target + "\t" + 100 + "\n";
	}

	public static void Uniprot2StringMap(String specie, String specieNo, String poslistfile, String backlistfile,
			String postxtFilePath, String backtxtFilePath, String outputFile) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outputFile)));
			BufferedReader br = new BufferedReader(new FileReader(new File(poslistfile)));
			LinkedHashMap<String, String> list = new LinkedHashMap<String, String>();
			String tmpstr = "";
			while ((tmpstr = br.readLine()) != null) {
				tmpstr = tmpstr.trim();
				if (tmpstr == "")
					continue;
				list.put(tmpstr, "");
				bw.write(getoneFileUniprot2String(specie, specieNo, tmpstr, postxtFilePath));
			}
			br = new BufferedReader(new FileReader(new File(backlistfile)));
			while ((tmpstr = br.readLine()) != null) {
				tmpstr = tmpstr.trim();
				if (tmpstr == "" || list.containsKey(tmpstr))
					continue;
				bw.write(getoneFileUniprot2String(specie, specieNo, tmpstr, backtxtFilePath));
			}
			bw.flush();
			bw.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	/**
	 * 将uniprot下载的fasta文件读入到程序中并输出UniprotID 注：需要调整java heap大小
	 * 
	 * @param infile
	 *            读入的文件
	 * @param outfile
	 *            输出的文件
	 * @throws Exception
	 */
	public static void getAcc(String infile, String outfile) throws Exception {
		LinkedHashMap<String, ProteinSequence> tmpmap = FastaReaderHelper.readFastaProteinSequence(new File(infile));
		String key[] = tmpmap.keySet().toArray(new String[0]);
		FileWriter fw = new FileWriter(new File(outfile));
		for (int i = 0; i < key.length; i++) {
			fw.write(key[i] + "\n");
		}
		fw.flush();
		fw.close();
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String spename = "HUMAN";
		FunctionalFeaturesUtil_1 getnm=new FunctionalFeaturesUtil_1();
		getnm.test(spename,getnm,0.05);
	}

	public void test(String spename,FunctionalFeaturesUtil_1 getnm,double mypvalue) {
		String root = "E:\\WM\\GOandKEGG\\data\\";
		System.out.println(root);
		String root2 = "E:\\WM\\tts9\\workspace\\DAP\\downloads\\2018_08_02\\classfyResult\\"+spename; // change
																	// here,
																	// 指定对应物种的TXT文件的目录
		try {
			GetPhos.creatFasta("C-linked (Man) tryptophan.",spename,root);//获取Phosphorylation文件
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		FunctionalFeatures ff2 = new FunctionalFeatures(root2 + "\\backgroundList_"+spename+".txt",
				root + "Phosphorylation_"+spename+".txt", mypvalue); // 第一个参数是所有的人类蛋白质的ACC,
															// 第二个参数是所有人类的Phosphorylation的ACC
															// 第三个参数是pvalue，这里pvalue=0.05,可以运行用户设置为其他值

		TXTExtractHelper helper2 = new TXTExtractHelper(ff2.backgroundProtMap, root2);
		helper2.extractandOut(root + "output/");

		ff2.addGo(root + "output/GoCTerms.txt", root + "output/GoFTerms.txt", root + "output/GoPTerms.txt");
		ff2.addKEGG(root + "KEGG/hsa.txt", root + "output/KEGGMap.txt");
		ff2.docalculate();
		getnm.getn=ff2.mappingn;
		getnm.getm=ff2.mappingm;
		System.out.println(getnm.getn+"-------"+getnm.getm);
		System.out.println("fine!");
	}
}
