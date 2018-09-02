package team.ptm.model.service;

import team.ptm.model.entity.*;
import java.util.List;
import org.apache.poi.hssf.usermodel.HSSFCell;
import org.apache.poi.hssf.usermodel.HSSFRow;
import org.apache.poi.hssf.usermodel.HSSFSheet;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.poifs.filesystem.POIFSFileSystem;

import org.apache.tools.zip.ZipEntry;
import org.apache.tools.zip.ZipOutputStream;

import net.sf.json.JSONArray;

/*import java.awt.List;*/
import java.io.*;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;

// 从TXT文件中获取注释的 N-linked (GlcNAc...) asparagine. 糖基化位点信息
/*
 * inputFilePath为需要获取信息的某个物种的uniprot IDs文件路径，文件中每行一个uniprot ID
 * annotationFileFolder为注释文件txt目录
 * finalResult为结果文件路径
 */
public class FindResult {

	private static String inputFilePath = "E:\\WM\\tts9\\workspace\\DAP\\downloads\\2018_08_02\\classfyResult";

	public static File getResultfiles() {
		return resultfiles;
	}

	public static void setResultfiles(File resultfiles) {
		FindResult.resultfiles = resultfiles;
	}

	private static File resultfiles;

	// 创建txt文件
	public static String getAnnotation(JSONArray specielist, JSONArray RqptmList, String serverDir) throws IOException {

		// 将结果文件上传到服务器上
		File uploadfile = new File(serverDir);
		if (!uploadfile.exists()) {
			boolean isCreated = uploadfile.mkdir();
			if (!isCreated) {
				// 目标上传目录创建失败,可做其他处理,例如抛出自定义异常等,一般应该不会出现这种情况。
				return null;
			}
		}

		// 结果文件以请求日期命名
		Date date = new Date();
		SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd HH-mm-ss");
		String resultstr = "result" + sdf.format(date);
		resultfiles = new File(serverDir, resultstr);
		resultfiles.mkdir();

		try {
			for (int m = 0; m < specielist.size(); m++) {
				System.out.println("dir is " + inputFilePath + File.separator + specielist.get(m).toString().trim());
				// 提前处理好的物种文件夹
				File fileDir = new File(inputFilePath + File.separator + specielist.get(m).toString().trim());
				System.out.println("result is " + resultfiles.getPath());

				if (fileDir.exists()) {
					/*
					 * File resultfile = new File(resultfiles.getPath(),
					 * specielist.get(m).toString().trim()); resultfile.mkdir();
					 */
					File finalResult = new File(resultfiles.getPath(), specielist.get(m).toString().trim() + ".txt");
					// 这条语句存在问题，返回的物种文件名不应该是标识码，而是物种名才对
					File[] files = fileDir.listFiles();

					OutputStreamWriter write = new OutputStreamWriter(new FileOutputStream(finalResult), "utf-8");
					BufferedWriter bw1 = new BufferedWriter(write);
					bw1.write("UniProtID\tPosition\tPTMType\tResidue\tPeptide\n");

					// 用list暂存结果
					List<Ptmfilter> tmpList = new ArrayList<Ptmfilter>();

					for (File tempfile : files) {
						String tmpstr = "";
						int flag = 0;
						int flag0 = 0;
						String[] str = new String[20];

						// 这三个变量用于处理sequence序列
						String[] tmpStr = new String[6];
						List<String> seqList = new ArrayList<String>();
						List<String> tseqList = new ArrayList<String>();
						BufferedReader br = new BufferedReader(new FileReader(tempfile));
						String Acc = getACC(tempfile, ".txt");
						String ptm = "";
						// 逐个读取注释文件
						while ((tmpstr = br.readLine()) != null) {
							// 去掉前后空格
							tmpstr = tmpstr.trim();

							// 第一行的糖基化信息

							if (flag == 0 && tmpstr.startsWith("FT   CARBOHYD")) {
								
								for (int rpos = 0; rpos < RqptmList.size(); rpos++) {
									String tstr = RqptmList.get(rpos).toString();
									if (tmpstr.endsWith(tstr)) {
										System.out.println(tmpstr);
										str = tmpstr.trim().split("\\s+");
										ptm = tstr;
										flag = 1;
									}
								}

							}
							// 第二行的ECO信息
							else if (flag == 1 && tmpstr.contains("ECO:0000269")) {

								System.out.println(Acc + "_" + str[2] + "_" + ptm);
								Ptmfilter tmpPtm = new Ptmfilter();
								tmpPtm.setUniprotID(Acc);
								tmpPtm.setPosition(Integer.parseInt(str[2]));
								tmpPtm.setPtmType(ptm);

								tmpPtm.setFlag(false);
								tmpList.add(tmpPtm);
								flag = 0;
							}
							// 为了得到位点信息，此处需要做标记
							if (tmpstr.startsWith("SQ")) {
								flag0 = 1;
								continue;
							}
							if (tmpstr.startsWith("//")) {
								flag0 = 0;
							}
							// 此处开始处理蛋白质的序列信息，
							// 根据已经获得的位点，将氨基酸以及所在序列段（8位）获取到

							if (flag0 == 1) {
								// list中最后一个ptmfilter的Flag
								// 为true，说明已经将位点等信息补全，可以跳出循环，继续打开下一个文件。
								if (tmpList.size() == 0) {
									flag0 = 0;
									break;
								}

								else if (tmpList.get(tmpList.size() - 1).getFlag()) {
									flag0 = 0;
									break;
								}
								// 将sequence遍历，字符串->去掉前后空格->字符串数组->转换成字符数组
								else {
									tmpStr = tmpstr.trim().split("\\s+");
									tseqList = Arrays.asList(tmpStr);
									seqList.addAll(tseqList);
								}

							}

						} // 一个文件读完

						// 回溯tmpList
						if (tmpList.size() == 0) {
							flag0 = 0;
							continue;
						}

						for (int p = tmpList.size() - 1; p >= 0; p--) {
							if (tmpList.get(p).getFlag()) {
								System.out.println("flag0 is " + flag0);
								break;
							} else {
								int pos = tmpList.get(p).getPosition();
								int listpos = (int) pos / 10;
								int abspos = pos % 10;
								String peptide = "";
								String pre = "";
								String end = "";
								if (abspos == 0) {
									tmpList.get(p).setResidue(seqList.get(listpos - 1).charAt(9));
									
									peptide = seqList.get(listpos).substring(3);
									if (seqList.get(listpos+1) != null){
										end = seqList.get(listpos+1);
										if(end.length()>6)
											end = end.substring(0, 6);
									}
									peptide = peptide+end;
									tmpList.get(p).setPeptide(peptide);
								} else {
									tmpList.get(p).setResidue(seqList.get(listpos).charAt(abspos - 1));
									
									if (abspos - 1 < 4) {
										int b0 = 10 - (6 - (abspos - 1));
										int e1 = abspos + 6;
										pre = seqList.get(listpos - 1).substring(b0);
										end = seqList.get(listpos);
										if (end.length() > (6 + abspos))
											end = end.substring(0, e1);
										peptide = pre + end;
										tmpList.get(p).setPeptide(peptide);
									}

									else if (abspos - 1 > 5) {
										int b1 = abspos - 1 - 6;
										int e2 = 6 - (9 - (abspos - 1));

										if (seqList.get(listpos + 1) != null) {
											pre = seqList.get(listpos).substring(b1);
											end = seqList.get(listpos + 1);
											if (end.length() > e2)
												end = end.substring(0, e2);
										} else {
											pre = seqList.get(listpos).substring(b1);
											end = "";
										}

										peptide = pre + end;
										tmpList.get(p).setPeptide(peptide);
									} else {
										int b0 = 10 - (6 - (abspos - 1));
										// int e2 =
										// 10+(6-(abspos-1))+(6-(9-(abspos-1)));
										int e2 = 13;
										peptide = seqList.get(listpos - 1) + seqList.get(listpos);
										if (seqList.get(listpos + 1) != null) {
											peptide += seqList.get(listpos + 1);
											peptide = peptide.substring(b0);
											if (peptide.length() > 13)
												peptide = peptide.substring(0, e2);
										}
										tmpList.get(p).setPeptide(peptide);
									}

								}

							}
							tmpList.get(p).setFlag(true);
							;
						}
						br.close();
					}

					// 将物种的所选的ptmlist的过滤信息存储在“物种名.txt”文件中。
					for (int t = 0; t < tmpList.size(); t++) {
						System.out.println(tmpList.get(t).getResidue());
						bw1.write(tmpList.get(t).getUniprotID() + "\t" + tmpList.get(t).getPosition() + "\t"
								+ tmpList.get(t).getPtmType() + "\t" + tmpList.get(t).getResidue() + "\t"
								+ tmpList.get(t).getPeptide()  + "\n");
						// 此处可以同时将物种在ptmlist的过滤后的蛋白质名字输入到文件“ptm_物种名.txt”中
						// 问题：去冗余和富集性分析输入的文件是一个ptm对应的物种蛋白质信息，还是所选择的所有ptm对应的物种蛋白质信息
					}

					bw1.close();

				} else {
					System.out.println(specielist.get(m).toString().trim() + " is not existed!");
				}
			}

			// 测试
			for (int m = 0; m < RqptmList.size(); m++) {
				System.out.println("RqptmList is " + RqptmList.get(m));
			}

		} catch (Exception e) {
			e.printStackTrace();
		}
		return resultfiles.getPath();
	}

	// 创建excel文件
	public static File txtToExcel(File ifile, String newDirpath) {
		String ofilename = getACC(ifile, ".txt");
		File ofile = new File(newDirpath + "//" + ofilename + ".xls");
		try {
			/* File test = new File("E:\\data\\test.txt"); */
			BufferedReader br = new BufferedReader(new FileReader(ifile));
			HSSFWorkbook hwb = new HSSFWorkbook();
			HSSFSheet sheet = hwb.createSheet("testexcel");
			HSSFRow firstrow = sheet.createRow(0);
			firstrow.setHeightInPoints(30);
			HSSFCell[] firstcell = null;
			String str = null;

			int i = 0;
			OutputStream out = new FileOutputStream(ofile);
			while ((str = br.readLine()) != null) {
				firstrow = sheet.createRow(i++);
				String[] list = str.split("\t");
				firstcell = new HSSFCell[list.length];
				for (int k = 0; k < list.length; k++) {
					System.out.println(list[k]);
					firstcell[k] = firstrow.createCell(k);
					firstcell[k].setCellValue(list[k]);
				}
				// firstrow = sheet.createRow(i);
				// firstcell = new HSSFCell[1];
				// firstcell[i] = firstrow.createCell(0);
				// firstcell[i].setCellValue(str);

			}

			hwb.write(out);
			out.close();
			br.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return ofile;
	}

	// 创建csv格式
	public static File excelToCSV(File ifile, String newDirpath) {
		String ofilename = getACC(ifile, ".xls");
		File ofile = new File(newDirpath + "//" + ofilename + ".csv");
		String buffer = "";

		try {
			// 1.读取Excel的对象
			POIFSFileSystem poifsFileSystem = new POIFSFileSystem(new FileInputStream(ifile));
			// 2.Excel工作薄对象
			HSSFWorkbook hssfWorkbook = new HSSFWorkbook(poifsFileSystem);
			// 3.Excel工作表对象
			HSSFSheet hssfSheet = hssfWorkbook.getSheetAt(0);
			// 总行数
			int rowLength = hssfSheet.getLastRowNum() + 1;
			// 4.得到Excel工作表的行
			HSSFRow hssfRow = hssfSheet.getRow(0);
			// 总列数
			int colLength = hssfRow.getLastCellNum();
			// 得到Excel指定单元格中的内容
			// HSSFCell hssfCell = hssfRow.getCell(0);
			// //得到单元格样式
			// CellStyle cellStyle = hssfCell.getCellStyle();

			for (int i = 0; i < rowLength; i++) {
				// 获取Excel工作表的行
				HSSFRow hssfRow1 = hssfSheet.getRow(i);
				for (int j = 0; j < colLength; j++) {
					// 获取指定单元格
					HSSFCell hssfCell1 = hssfRow1.getCell(j);

					// Excel数据Cell有不同的类型，当我们试图从一个数字类型的Cell读取出一个字符串时就有可能报异常：
					// Cannot get a STRING value from a NUMERIC cell
					// 将所有的需要读的Cell表格设置为String格式
					if (hssfCell1 != null) {
						System.out.print(hssfCell1.getStringCellValue() + "\t");
						buffer += hssfCell1.getStringCellValue().replaceAll("\n", " ") + ",";
					}
				}
				buffer = buffer.substring(0, buffer.lastIndexOf(",")).toString();
				buffer += "\n";
			}
		} catch (Exception e) {
			e.printStackTrace();
		}

		// write the string into the file

		try {
			if (!ofile.exists())
				ofile.createNewFile();
			BufferedWriter writer = new BufferedWriter(new FileWriter(ofile));
			writer.write(buffer);
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		if (ofile.exists()) {
			System.out.println(ofile.getAbsolutePath());
		} else
			System.out.println("ofile is not existed");
		//deleteFile(ifile);
		return ofile;
	}

	private static String getACC(File temp, String end) {
		String acc = temp.getName();
		acc = acc.substring(0, acc.indexOf(end));
		// System.out.println("去掉后缀的文件名是："+acc);
		return acc;
	}

	public static boolean deleteFile(File file) {

		// 如果文件路径所对应的文件存在，并且是一个文件，则直接删除
		if (file.exists() && file.isFile()) {
			if (file.delete()) {
				System.out.println("删除单个文件" + file.getName() + "成功！");
				return true;
			} else {
				System.out.println("删除单个文件" + file.getName() + "失败！");
				return false;
			}
		} else {
			System.out.println("删除单个文件失败：" + file.getName() + "不存在！");
			return false;
		}
	}

	public static File compressFile(File zipFile, String desFile) {
		OutputStream outputStream = null;
		ZipOutputStream zipOutputStream = null;
		File desZipfile = new File(desFile);
		try {
			outputStream = new FileOutputStream(desZipfile);
			zipOutputStream = new ZipOutputStream(outputStream);
			zipOutputStream.setEncoding("utf-8");
			for (File file : zipFile.listFiles()) {
				startZip(zipOutputStream, "", file);
			}
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			try {
				if (zipOutputStream != null) {
					zipOutputStream.close();
				}
				if (outputStream != null) {
					outputStream.close();
				}
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		return desZipfile;
	}

	// 以防文件里面套文件
	private static void startZip(ZipOutputStream zos, String s, File file) {
		try {
			System.out.println("当前压缩文件： " + file.getName());
			if (file.isDirectory()) {
				System.out.println("正在压缩文件夹");
				String newpath1 = s + file.getName() + "/";
				ZipEntry ze1 = new ZipEntry(newpath1);
				zos.putNextEntry(ze1);
				zos.closeEntry();
				File[] files = file.listFiles();
				for (File f : files) {
					if (f.isDirectory()) {
						/*
						 * String newpath2 = newpath1+f.getName()+"/"; ZipEntry
						 * ze2 = new ZipEntry(newpath2);
						 */
						zos.putNextEntry(ze1);
						zos.closeEntry();
						startZip(zos, newpath1, f);
					} else {
						zipFile(zos, newpath1, f);
					}
				}
			} else {
				zipFile(zos, s, file);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private static void zipFile(ZipOutputStream zos, String s, File f) {
		System.out.println("正在压缩文件：" + f.getName());
		InputStream inputStream = null;
		BufferedInputStream bin = null;
		try {
			zos.putNextEntry(new ZipEntry(s + f.getName()));
			inputStream = new FileInputStream(f);
			bin = new BufferedInputStream(inputStream);
			byte[] buffer = new byte[bin.available()];
			Integer length = null;
			while ((length = bin.read(buffer)) != -1) {
				zos.write(buffer, 0, length);
			}
			zos.closeEntry();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static File getRsFile(String format, String ifcompress, String tempRsPath, String uploadDir) {
		// TODO Auto-generated method stub
		/*
		 * FASTA</option> <option value="txt" selected>TXT</option> <option
		 * value="excel">EXCEL</option> <option value="esv">ESV</op
		 */
		String rsurl = "";
		File rsDir = null;
		File rsZip = null;

		switch (format) {
		case "txt":
			if (ifcompress.equals("yes")) {
				// 压缩
				// 结果文件以请求日期命名
				Date date = new Date();
				SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd HH-mm-ss");
				String resultstr = "result" + sdf.format(date);

				rsDir = new File(uploadDir, resultstr);
				rsDir.mkdir();
				File tempRsfile = new File(tempRsPath);
				File RsDir = txtRenameDir(tempRsfile.getPath(), rsDir.getPath());
				// 压缩
				// 定义打包文件zip的存放位置
				String outputpath = rsDir.getAbsolutePath() + ".zip";
				// 使用ZIP工具类来压缩zipFile文件夹
				rsZip = FindResult.compressFile(RsDir, outputpath);
				System.out.println("压缩成功：" + rsZip.getAbsolutePath());
				rsurl = rsZip.getAbsolutePath();
			} else {
				// 直接显示文件内容
			}
			break;

		case "excel":
			if (ifcompress.equals("yes")) {
				// 压缩
				// 结果文件以请求日期命名
				Date date = new Date();
				SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd HH-mm-ss");
				String resultstr = "result" + sdf.format(date);

				rsDir = new File(uploadDir, resultstr);
				rsDir.mkdir();
				File tempRsfile = new File(tempRsPath);
				File RsDir = excelRenameDir(tempRsfile.getPath(), rsDir.getPath());
				// 压缩
				// 定义打包文件zip的存放位置
				String outputpath = rsDir.getAbsolutePath() + ".zip";
				// 使用ZIP工具类来压缩zipFile文件夹
				rsZip = FindResult.compressFile(RsDir, outputpath);
				System.out.println("压缩成功：" + rsZip.getAbsolutePath());
				rsurl = rsZip.getAbsolutePath();
			} else {
				// 直接显示文件内容
			}
			break;
		case "esv":
			if (ifcompress.equals("yes")) {
				// 压缩
				// 结果文件以请求日期命名
				Date date = new Date();
				SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd HH-mm-ss");
				String resultstr = "result" + sdf.format(date);

				rsDir = new File(uploadDir, resultstr);
				rsDir.mkdir();
				File tempRsfile = new File(tempRsPath);
				File RsDir = esvRenameDir(tempRsfile.getPath(), rsDir.getPath());
				// 压缩
				// 定义打包文件zip的存放位置
				String outputpath = rsDir.getAbsolutePath() + ".zip";
				// 使用ZIP工具类来压缩zipFile文件夹
				rsZip = FindResult.compressFile(RsDir, outputpath);
				System.out.println("压缩成功：" + rsZip.getAbsolutePath());
				rsurl = rsZip.getAbsolutePath();

			} else {
				// 直接显示文件内容
			}
			break;
		case "fasta":
			if (ifcompress.equals("yes")) {
				// 需要进一步处理

			} else {
				// 直接显示文件内容
			}
			break;

		default:
			break;
		}

		return rsZip;

	}

	private static File esvRenameDir(String oldDirpath, String newDirpath) {
		// TODO Auto-generated method stub
		// TODO Auto-generated method stub
		File fl = new File(oldDirpath); // 这里写上发替换的文件夹路径,注意使用双斜杠
		String[] files = fl.list();
		File newDir = new File(newDirpath);
		newDir.mkdir();
		File f = null;
		File excelF = null;

		for (String file : files) {
			f = new File(fl, file);// 注意,这里一定要写成File(fl,file)如果写成File(file)是行不通的,一定要全路径
			excelF = txtToExcel(f, newDirpath);
			excelToCSV(excelF, newDirpath);
		}
		return newDir;
	}

	private static File excelRenameDir(String oldDirpath, String newDirpath) {
		// TODO Auto-generated method stub
		// TODO Auto-generated method stub
		File oldDir = new File(oldDirpath); // 这里写上发替换的文件夹路径,注意使用双斜杠
		String[] files = oldDir.list();
		File newDir = new File(newDirpath);
		newDir.mkdir();

		File f = null;
		String filename = "";
		for (String file : files) {
			f = new File(oldDir, file);// 注意,这里一定要写成File(fl,file)如果写成File(file)是行不通的,一定要全路径

			txtToExcel(f, newDirpath);
			// System.out.println(filename);
			/*
			 * the string 要替换掉的内容 is the content in your own file string with
			 * the name 替换成的内容, here you should change the string into what you
			 * have.
			 */

		}
		return newDir;
	}

	private static File txtRenameDir(String oldDirpath, String newDirpath) {
		// TODO Auto-generated method stub
		File fl = new File(oldDirpath); // 这里写上发替换的文件夹路径,注意使用双斜杠
		String[] files = fl.list();
		File newDir = new File(newDirpath);
		newDir.mkdir();
		File f = null;
		String filename = "";
		for (String file : files) {
			f = new File(fl, file);// 注意,这里一定要写成File(fl,file)如果写成File(file)是行不通的,一定要全路径
			filename = f.getName();
			// System.out.println(filename);
			/*
			 * the string 要替换掉的内容 is the content in your own file string with
			 * the name 替换成的内容, here you should change the string into what you
			 * have.
			 */
			File newFile = new File(newDirpath + File.separator + filename);
			f.renameTo(newFile);
		}
		return newDir;
	}
}