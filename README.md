1. 獲取Strandard cell的Matel layers的多邊形位置資訊

準備前 :

請在當前資料夾放入以下資料
get_polys.il(自行準備但資料夾有範例，需將檔案內libName修改為當前Library命名、outStream修改為輸出的檔案位置與名稱)

執行(Virtuoso) :
load("get_polys.il")

輸出 :
Matel layers的多邊形位置資訊(預設為.csv檔)

2. 根據cell相對位置與內部Matel layers分布判斷是否可交換以最大化繞線資源

準備前 :

請在當前資料夾放入以下資料
(1) TEST_EVEN.csv、TEST_ODD.csv...(各種Library內Strandard cell的Matel layers位置資訊，由前面get_polys.il所生成)
(2) PATH.csv...(自行準備但資料夾有範例，根據產生的Library位置資訊檔(前面)執行時的位置做更動)
(注 : 此程式包含讓Cell於APR時，使pin始終On track的優化，順序依Strandard cell的matel範圍由大至小排列，以其中TIELOx1為例:)
****************
cell TIELOx1
-1 0
-2 1
0 0
1 1
****************
為了使pin始終On track，在此範例會區分EVEN、ODD兩種Library模式，根據Cell在APR時的位置做判斷修改
在此-2為EEVN、-1為ODD，他們後面必須有至少一個Library，在此範例為-1對上0(執行時第一個位置)、-2對上1(執行時第二個位置)，以此類推根據Strandard cell自身的matel範圍由大至小排列
之後每個Library也要對應有它的下一步，以0為例，其下一步先為自身0，此例EVEN/ODD各只有一個，若還有能接續在後面(ex: 0 0 2 3...)
(3) output.def...(自行準備但資料夾有範例，於ICC2使用write_def所生成)

編譯 :
g++ -o CHANGE CHANGE.cpp

執行 :
./CHANGE TEST_EVEN.csv TEST_ODD.csv... PATH.csv output.def 
(注 : ... 位置表示可根據需要放置多個library，PATH.csv需完整配合調整，放置越多可做的調整越精細)

輸出 :
output.tcl (內部含有需做調整的cells的tcl指令，於APR(ICC2)的CTS後Routing前下此指令完成更改，以最大化繞線資源)
