; Specify library name and output file
let((libName outStream cellList cellName dbcv M1 bboxList)
    libName = "INVx1" ; Replace with your library name
    outStream = outfile("/home/s1126002/TEST/APR/CH/INVx1.csv")

    ; Get all cells in the library
    cellList = ddGetObj(libName)~>cells

    ; Iterate through each cell
    foreach(cell cellList
        cellName = cell~>name
        printf("Processing cell: %s\n" cellName)
		fprintf(outStream "cell %s " cellName)

        ; Open the cell view in read mode
        dbcv = dbOpenCellViewByType(libName cellName "layout" "" "r")
        
        ; Ensure the cell view exists
        when(dbcv
			; Query BOUNDARY shapes
            BOUNDARY = dbShapeQuery(dbcv list("BOUNDARY" "drawing") (dbcv~>bBox))
			
			; Get bounding boxes for BOUNDARY
            foreach(shape BOUNDARY
                let((bbox)
                    bbox = shape~>bBox
                    ; Write BOUNDARY size to the CSV file
                    fprintf(outStream "%f %f\n"
                            (car (cadr bbox)) ; MaxX
                            (cadr (cadr bbox)) ; MaxY
                    )
                )
            )
			
            ; Query M1 shapes
            M1 = dbShapeQuery(dbcv list("M1" "drawing") (dbcv~>bBox))
			; Merge M1
			M1merge = leMergeShapes(M1)

            ; Get bounding boxes for M1
            foreach(shape M1merge
                let((pts bbox)
                    pts = shape~>points
					bbox = shape~>bBox
                    ; Write M1 coordinates to the CSV file
					; if pts is exist, use pts to output
					foreach(pt pts   
						fprintf(outStream "( %f %f )," (car pt) (cadr pt))
					)
					; if pts is nil, it is Rect use bbox to output
					when(!pts
						fprintf(outStream "( %f %f ),( %f %f ),( %f %f ),( %f %f ),"
							(car (cadr bbox))(cadr (cadr bbox)); MaxX MaxY
							(car (car bbox))(cadr (cadr bbox)); MinX MaxY
							(car (car bbox))(cadr (car bbox)); MinX MinY
							(car (cadr bbox))(cadr (car bbox)); MaxX MinY
						)
                    )
                )
				fprintf(outStream "\n")
            )
            ; Close the cell view
            dbClose(dbcv)
        )
    )

    ; Close the output file
    close(outStream)

    ; Completion message
    printf("Export complete")
)
