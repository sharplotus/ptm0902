package team.ptm.model.service;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import org.apache.poi.hssf.usermodel.HSSFCell;
import org.apache.poi.hssf.usermodel.HSSFRow;
import org.apache.poi.hssf.usermodel.HSSFSheet;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import team.ptm.model.entity.mapping;
import team.ptm.model.entity.mappinggo;

/**
 * @author ivy 2018.8.27
 * 
 *         调用命令行执行jar文件获得图片
 *
 */
public class EnrichPic {

	// 获取图片需要csv文件和R文件，csv需要通过txt到xls的转换，R文件需要文件内容的读写，最后用命令行执行R文件
	public static int rowkey = 0;

	public static void main(String[] args) {
		// 生成mapping前使用的文件
		String root = "E:/WM/GOandKEGG/data/GO/term/";
		String spename = "HUMAN";

		FunctionalFeaturesUtil_1 getnm = new FunctionalFeaturesUtil_1();
		getnm.test(spename, getnm, 0.05);

		// 生成mapping文件
		System.out.println(root);
		try {
			FunctionMapping.GOMapping(root + "GoC.txt", root + "mapping/" + spename + "mapping_GoC.txt");
			FunctionMapping.GOMapping(root + "GoP.txt", root + "mapping/" + spename + "mapping_GoP.txt");
			FunctionMapping.GOMapping(root + "GoF.txt", root + "mapping/" + spename + "mapping_GoF.txt");
			FunctionMapping.KEGG_Mapping(getnm, spename, root + "KEGG.txt",root + "mapping/" + spename + "mapping_KEGG.txt");
		} catch (Exception e) {
			e.printStackTrace();
		}
		// 将mapping.txt转化成xls文件
		// 在xls文件中处理数据
		// 将xls文件转化成csv文件
		filechange(root + "mapping/" + spename + "mapping_KEGG.txt", root + "mapping/");
		// 生成需要的R文件
		createR(root + "mapping/" + spename + "mapping_KEGG.R", root + "mapping/" + spename + "mapping_KEGG.csv",
				spename);
		System.out.println(root + "mapping/" + spename + "mapping_KEGG.csv");
		// 寻找图片
		gainPic("\"" + root + "mapping/" + spename + "mapping_KEGG.R\"");

		// 以下为GO文件的成图过程

		File agc_akt = new File(root + "mapping/" + spename + "mapping_agcakt.xls");

		List<File> goList = new ArrayList<File>();
		goList.add(FindResult.txtToExcel(new File(root + "mapping/" + spename + "mapping_GoP.txt"), root + "mapping/"));
		goList.add(FindResult.txtToExcel(new File(root + "mapping/" + spename + "mapping_GoC.txt"), root + "mapping/"));
		goList.add(FindResult.txtToExcel(new File(root + "mapping/" + spename + "mapping_GoF.txt"), root + "mapping/"));
		try {
			readGOExcel(goList, agc_akt);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		FindResult.excelToCSV(agc_akt, root + "mapping/");
	}

	/**
	 * csv需要通过txt到xls的转换
	 * 
	 * @param txtfiel
	 *            输入文件，原始的txt类型的mapping文件
	 * @param mappath
	 *            存放cvs文件的路径
	 */
	public static void filechange(String txtfielname, String mappath) {

		File txtfiel = new File(txtfielname);
		// xlsfile转化后的xls文件
		File xlsfile = FindResult.txtToExcel(txtfiel, mappath);// 从txt文件转化成xls文件
		try {
			readExcel(xlsfile);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} // 使用中间生成的xls文件处理数据
		FindResult.excelToCSV(xlsfile, mappath);// 从xls文件转化成scv文件
	}

	/**
	 * 对xls文件进行排序
	 * 
	 * @param filename
	 * @return
	 * @throws Exception
	 */
	public static void readExcel(File filename) throws Exception {
		ArrayList<mapping> map = new ArrayList<mapping>();

		try {
			InputStream in = new FileInputStream(filename);
			HSSFWorkbook workbook = new HSSFWorkbook(in);
			HSSFSheet sheet = workbook.getSheetAt(0); // sheet 从0开始
			int rowNum = sheet.getLastRowNum() + 1; // 取得最后一行的行号

			for (int i = 1; i < rowNum; i++) { // 行循环开始，从第2行开始，索引为1
				mapping temp = new mapping();
				HSSFRow row = sheet.getRow(i); // 行
				// System.out.println(row.getCell((short) 1));//这是第2列，索引为1
				if (row == null)
					continue; // 中间如果有空行，则退出
				// 将mapping的每一行加入到
				temp.setID(row.getCell((short) 0).toString());
				temp.setCount(Integer.parseInt(row.getCell((short) 1).toString()));
				temp.setM(Integer.parseInt(row.getCell((short) 2).toString()));
				temp.setPvalue(Double.parseDouble(row.getCell((short) 3).toString()));
				temp.setTerm(row.getCell((short) 4).toString());
				temp.setN_dbz(Integer.parseInt(row.getCell((short) 5).toString()));
				temp.setM_plus(Integer.parseInt(row.getCell((short) 6).toString()));
				map.add(temp);
			}
			in.close();
			SortByPvalue(map);

			// System.out.println("----------------------------");
			int count = 0;
			HSSFWorkbook workbook2 = new HSSFWorkbook();
			HSSFSheet sheet2 = workbook2.createSheet();
			sheet2 = workbook2.getSheetAt(0);

			HSSFRow row2 = sheet2.createRow(0);
			row2.createCell(0).setCellValue("ID");
			row2.createCell(1).setCellValue("count");
			row2.createCell(2).setCellValue("M");
			row2.createCell(3).setCellValue("Pvalue");
			row2.createCell(4).setCellValue("Term");
			row2.createCell(5).setCellValue("n");
			row2.createCell(6).setCellValue("m");

			for (int i = 0; i < map.size(); i++) { // 行循环开始，从第2行开始，索引为1
				if (count >= 20) {// 当数据量大于20时舍弃
					break;
				}
				row2 = sheet2.createRow(i + 1);
				row2.createCell(0).setCellValue(map.get(i).getID());
				row2.createCell(1).setCellValue(String.valueOf(map.get(i).getCount()));
				row2.createCell(2).setCellValue(String.valueOf(map.get(i).getM()));
				row2.createCell(3).setCellValue(String.valueOf(map.get(i).getPvalue()));
				row2.createCell(4).setCellValue(map.get(i).getTerm());
				row2.createCell(5).setCellValue(String.valueOf(map.get(i).getN_dbz()));
				row2.createCell(6).setCellValue(String.valueOf(map.get(i).getM_plus()));

				count++;
			}
			FileOutputStream out;
			out = new FileOutputStream(filename);
			workbook2.write(out);// 保存Excel文件
			out.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println("mapping中数据处理完成！");
	}

	/**
	 * 对xls文件进行排序
	 * 
	 * @param golist
	 *            一个xls文件的list包含三个go文件
	 * @param agc_akt
	 *            要获得的处理后的包含三个go信息的xls文件
	 * @return
	 * @throws Exception
	 */
	public static File readGOExcel(List<File> golist, File agc_akt) throws Exception {
		// 写入的设置
		HSSFWorkbook workbook2 = new HSSFWorkbook();
		HSSFSheet sheet2;

		for (int f = 0; f < golist.size(); f++) {

			ArrayList<mappinggo> mapgo = new ArrayList<mappinggo>();
			File filename = golist.get(f);

			// 设置第一例
			String cate = null;
			if (filename.getName().contains("GoC")) {
				cate = "GOTERM_CC_FAT";
			} else if (filename.getName().contains("GoF")) {
				cate = "GOTERM_MF_FAT";
			} else if (filename.getName().contains("GoP")) {
				cate = "GOTERM_BP_FAT";
			}
			try {
				// 设置其他列
				InputStream in = new FileInputStream(filename);
				HSSFWorkbook workbook = new HSSFWorkbook(in);
				HSSFSheet sheet = workbook.getSheetAt(0);
				int rowNum = sheet.getLastRowNum() + 1;

				for (int i = 0; i < rowNum; i++) {
					mappinggo temp = new mappinggo();
					HSSFRow row = sheet.getRow(i);
					if (row == null)
						continue;
					temp.setCategory(cate);
					temp.setCount(String.valueOf(row.getCell((short) 1)));
					temp.setTerm(row.getCell((short) 4).toString().replaceAll(",", "  "));
					temp.setPValue(Double.parseDouble(row.getCell((short) 3).toString()));
					// 获得当前go文件的内容存到一个list里面
					mapgo.add(temp);
				}
				in.close();
				// 排好序
				GOSortByPvalue(mapgo);

				// 进行写
				// 寻找表页
				if (workbook2.getSheet("Sheet0") == null) {
					sheet2 = workbook2.createSheet("Sheet0");
					System.err.println("创建！！！！！！！！！！");
				} else {
					sheet2 = workbook2.getSheet("Sheet0");
					System.err.println("获取！！！！！！！！！！");
				}
				// 寻找行
				HSSFRow row2;
				int rowNum2 = sheet2.getLastRowNum();// 获得当前的最后一行
				if (rowNum2 == 0) {
					row2 = sheet2.createRow(0);
					row2.createCell(0).setCellValue("ID");
					row2.createCell(1).setCellValue("Count");
					row2.createCell(3).setCellValue("Term");
					row2.createCell(2).setCellValue("PValue");
				}
				rowNum2++;
				System.out.println("========" + rowNum2 + "==============================");
				int count = 0;// 设置前五
				// rowNum2 += rowkey;
				for (int i = rowNum2; i < rowNum2 + mapgo.size(); i++) { // 行循环开始，从第2行开始，索引为1
					if (count >= 5) {// 当数据量大于20时舍弃
						break;
					}
					System.out.println("========" + i);
					row2 = sheet2.createRow(i);
					row2.createCell(0).setCellValue(mapgo.get(i - rowNum2).getCategory());
					row2.createCell(1).setCellValue(mapgo.get(i - rowNum2).getCount());
					row2.createCell(2).setCellValue(String.valueOf(mapgo.get(i - rowNum2).getPValue()));
					row2.createCell(3).setCellValue(mapgo.get(i - rowNum2).getTerm());

					count++;
				}

			} catch (Exception e) {
				e.printStackTrace();
			}
			System.out.println("mappingdo*3中数据处理完成！" + filename.getName());
		}
		FileOutputStream out;
		out = new FileOutputStream(agc_akt);
		workbook2.write(out);// 保存Excel文件
		out.close();

		return agc_akt;
	}

	/**
	 * 给map文件排序:从小到大
	 * 
	 * @param map
	 */
	public static void SortByPvalue(ArrayList<mapping> map) {
		for (int i = 0; i < map.size() - 1; i++) {
			for (int j = 1; j < map.size() - i; j++) {
				mapping a = new mapping();
				if (map.get(j - 1).getPvalue() > map.get(j).getPvalue()) {
					a = map.get(j - 1);
					map.set((j - 1), map.get(j));
					map.set(j, a);
				}
			}
		}

	}

	/**
	 * 给mapgo文件排序:从小到大
	 * 
	 * @param mapgo
	 */
	public static void GOSortByPvalue(ArrayList<mappinggo> mapgo) {
		for (int i = 0; i < mapgo.size() - 1; i++) {
			for (int j = 1; j < mapgo.size() - i; j++) {
				mappinggo a = new mappinggo();
				if (mapgo.get(j - 1).getPValue() > mapgo.get(j).getPValue()) {
					a = mapgo.get(j - 1);
					mapgo.set((j - 1), mapgo.get(j));
					mapgo.set(j, a);
				}
			}
		}

	}

	/**
	 * 对R文件的操作,根据一个已有的固定模板创建一个当前需要的R文件用于画图
	 * 
	 * @param reafone
	 *            输入文件，一个包含R文件内容的txt文件
	 * @param writeone
	 *            输出文件，需要的新R文件
	 * @param cvsname
	 *            R画图使用到的cvs
	 */

	public static void createR(String writeone, String cvsname, String spename) {
		try {
			FileReader fr = new FileReader("E:/WM/templetR.txt");
			FileWriter fw = new FileWriter(writeone);
			PrintWriter pw = new PrintWriter(fw);
			System.out.println(cvsname);
			// cvsname.replaceAll("/", "//");
			pw.println("rose <- read.csv(\"" + cvsname + "\")");
			pw.println("setwd(\"E:/WM\")");
			pw.println("library(ggplot2)");
			pw.println("png(\"" + spename + "_KEGG.png\")");

			int read;
			read = fr.read();
			while (read != -1) {
				fw.write(read);
				read = fr.read();
			}
			fr.close();
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * 最后用命令行执行R文件
	 * 
	 * @param rroot
	 */
	public static void gainPic(String rroot) {
		Process pro;
		try {
			pro = Runtime.getRuntime().exec("\"C:/Program Files/R/R-3.5.1/bin/Rscript.exe\" " + rroot);
			System.out.println("已经生成图片");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
