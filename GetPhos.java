package team.ptm.model.service;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;

import net.sf.json.JSONArray;
import team.ptm.model.entity.Ptmfilter;

public class GetPhos {

	public static void main(String[] args) {
		try {
			creatFasta("C-linked (Man) tryptophan.","HUMAN","E:\\WM\\GOandKEGG\\data\\");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	private static String inputFilePath = "E:\\WM\\tts9\\workspace\\DAP\\downloads\\2018_08_02\\classfyResult";


	// 创建fasta文件
	public static List<String> creatFasta(String ptmStr,String spStr,String serverDir) throws IOException {
		List<String> rslist = new ArrayList<String>();
		// 将fasta文件暂存在本地服务器上
		//之后需要存储在服务器上：/tempPTM
		File uploadfile = new File(serverDir);
		if (!uploadfile.exists()) {
			boolean isCreated = uploadfile.mkdir();
			if (!isCreated) {
				// 目标上传目录创建失败,可做其他处理,例如抛出自定义异常等,一般应该不会出现这种情况。
				return null;
			}
		}
		String ptmStrname="Phosphorylation";
//		if(ptmStr.contains(" ")){
//			ptmStrname = ptmStr.replaceAll(" +", "_");
//			ptmStrname = ptmStrname.replace("(", "");
//			ptmStrname = ptmStrname.replace(")", "");
//
//		}
		//String filename1 = ptmStrname+"_"+spStr+".fasta";
		String filename2 = ptmStrname+"_"+spStr+".txt";
		//File file_fasta = new File(serverDir,filename1);//去冗余的输入文件
		File file_txt = new File(serverDir,filename2);//去冗余的输入文件
		//rslist.add(filename1);
		//rslist.add(file_fasta.getPath());
		rslist.add(ptmStr);
		try {
			
				System.out.println("dir is " + inputFilePath + File.separator + spStr.trim());
				// 提前处理好的物种文件夹
				File fileDir = new File(inputFilePath + File.separator + spStr.trim().trim());
				//System.out.println("result is " + file_fasta.getPath());

				if (fileDir.exists()) {
					
					// 这条语句存在问题，返回的物种文件名不应该是标识码，而是物种名才对
					File[] files = fileDir.listFiles();

					//OutputStreamWriter write1 = new OutputStreamWriter(new FileOutputStream(file_fasta), "utf-8");
					OutputStreamWriter write2 = new OutputStreamWriter(new FileOutputStream(file_txt), "utf-8");

					//BufferedWriter bw1 = new BufferedWriter(write1);
					BufferedWriter bw2 = new BufferedWriter(write2);

					for (File tempfile : files) {
						String tmpstr = "";
						int flag = 0;
						int flag0 = 0;
						int flag1 = 0;

						// 这三个变量用于处理sequence序列
						BufferedReader br = new BufferedReader(new FileReader(tempfile));
						String Acc = getACC(tempfile, ".txt");
						// 逐个读取注释文件
						while ((tmpstr = br.readLine()) != null) {
							// 去掉前后空格
							tmpstr = tmpstr.trim();

							if(flag1 == 0){
								// 第一行的糖基化信息
								if (flag == 0 && tmpstr.startsWith("FT   CARBOHYD")&&tmpstr.endsWith(ptmStr)) {
									System.out.println(tmpstr);
											flag = 1;
									}
								
								// 第二行的ECO信息
								else if (flag == 1 && tmpstr.contains("ECO:0000269")) {
									System.out.println(">sp|"+Acc+"| OS="+ spStr);
									//bw1.write(">sp|"+Acc+"|OS="+ spStr+"\n");
									bw2.write(Acc+"\n");
									flag = 0;
									flag1 = 1;
								}
							}
							else{
								// 为了得到序列信息，此处需要做标记
								if (tmpstr.startsWith("SQ")) {
									flag0 = 1;
									continue;
								}
								if (tmpstr.startsWith("//")) {
									flag0 = 0;
									flag1 = 0;
								}
								// 此处开始处理蛋白质的序列信息，
								if (flag0 == 1) {
									//bw1.write(tmpstr+"\n");
								}

							}
							
						}
					}
					//bw1.close();
					bw2.close();
				} else {
					System.out.println(spStr.trim() + " is not existed!");
				}
			

		} catch (Exception e) {
			e.printStackTrace();
		}
		return rslist;
	}
	
	private static String getACC(File temp, String end) {
		String acc = temp.getName();
		acc = acc.substring(0, acc.indexOf(end));
		// System.out.println("去掉后缀的文件名是："+acc);
		return acc;
	}
	
	

}
