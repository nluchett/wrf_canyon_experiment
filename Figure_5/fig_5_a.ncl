;   Example script to produce standard plots for a WRF quarter_ss run

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
;load "./WRFUserARW.ncl"

begin
;
; The WRF ARW input file.
; This needs to have a ".nc" appended, so just do it.

  a = addfile("/Volumes/Seagate Backup Plus Drive/WRF_output/Canyon_Experiment/flat_canyon_close/wrfout_d02_0001-01-01_06:00:30","r")

;  We generate plots, but what kind do we prefer?

;   type = "x11"
   type = "pdf"
;  type = "ps"
;  type = "ncgm"

  setvalues NhlGetWorkspaceObjectId()
    "wsMaximumSize" : 300000000
  end setvalues


  wks = gsn_open_wks(type,"flat_no_fire")
;  wks = gsn_open_wks(type,"plt_valley_vv_dbz_planar_D2")
;  gsn_define_colormap(wks,"WhViBlGrYeOrRe")
  gsn_define_colormap(wks,"BlWhRe")

; Set some basic resources
    res = True
    res@MainTitle     = "0~F34~0~F~ Slope"
    res@InitTime      = False
    res@Footer        = False
    res@FontHeightF   = 0.015
    res@MainTitlePos  = "Left"
    res@cnInfoLabelFontHeightF   = 0.02
    ;res@tiMainFontHeightF = 0.08
    ;res@gsnStringFontHeightF = 0.035   
    pltres = True
    pltres@FontHeightF   = 0.015
    plot_width = .6
    plot_height = .6

;  What times and how many time steps are in the data set?

  times  = wrf_user_list_times(a)  ; get times in the file
  ntimes = dimsizes(times)         ; number of times in the file


do it = 0,ntimes-10

    time = it
    res@TimeLabel                  = times(it)
    
    ;res@vpWidthF                   = plot_width
    ;res@vpHeightF                  = plot_height

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; First get the variables we will need

  p   = wrf_user_getvar(a, "pressure",time)    ; pressure
  th  = wrf_user_getvar(a,"th",time)    ; get temperature (C)
  u   = wrf_user_getvar(a,"ua",time)    ; ua is u averaged to mass points
  v   = wrf_user_getvar(a,"va",time)    ; va is v averaged to mass points
  w   = wrf_user_getvar(a,"W",time)    ; vertical velocity
  z   = wrf_user_getvar(a, "z",time)    ; grid point height
  ter = wrf_user_getvar(a,"HGT",time)   ; need terrain height sometimes
  dbz = wrf_user_getvar(a,"dbz",time)

;print(u(7,:,:))
printVarSummary(w)
printVarSummary(u)
printVarSummary(v)
 ; Vertical Velocity
      w_plane  = w(2,:,:)
      w_plane = w_plane
      opts_w   = res
      opts_w@FieldTitle            = "50 m Vertical Velocity"
      opts_w@UnitLabel             = "m s~S~-1~N~"
      ;opts_w@PlotLevelID           = 0.001*height + " km"
      opts_w@cnFillOn              = True
      opts_w@ContourParameters = (/ -5., 5., .1/)
      ;opts_w@gsnSpreadColorEnd     = -3
      opts_w@lbLabelFontHeightF = 0.020
      opts_w@lbTitleFontHeightF= .025
     ; opts_w@lbLabelStrings = (/"-5.0","-4.5","-4.0","-3.5","-3.0","-2.5","-2.0","-1.5","-1.0","-0.5","0","0.5","1.0","1.5","2.0","2.5","3.0","3.5","4.0","4.5","5.0"/)
      opts_w@tmYLMode          = "Explicit"     ; Set tick mark mode.
      opts_w@tmYLValues        = (/ 0, 100, 200, 300, 400 /)
      opts_w@tmYLLabels        = (/ -6000,  -3000, 0,  3000,  6000 /)
      opts_w@tmXBMode          = "Explicit"     ; Set tick mark mode.
      opts_w@tmXBValues        = (/ 0, 100, 200, 300, 400 /)
      opts_w@tmXBLabels        = (/ -6000,  -3000, 0,  3000,  6000 /)
      opts_w@tiYAxisString   = "Distance (m)"      ; Label for the Y axis
      opts_w@tiXAxisString   = "Distance (m)"      ; Label for the X axis
      
      contour_w = wrf_contour(a,wks, w_plane,opts_w)

  ; Winds
     u_plane  = u(0,:,:)
      v_plane  = v(0,:,:)
      opts_vct = res
      opts_vct@FieldTitle          = "Surface Winds"
      opts_vct@UnitLabel           = "m s~S~-1~N~"
      ;opts_vct@PlotLevelID         = height + " mb"
      opts_vct@NumVectors          = 20
      opts_vct@vcGlyphStyle        = "CurlVector"
      opts_vct@vcRefAnnoOn         = True
      opts_vct@vcRefMagnitudeF = 20
      opts_vct@vcRefAnnoFontHeightF = 0.023
      ;res@vcRefMagnitudeF = 4
      ;res@vcRefLengthF = .5
      ;res@vcMonoLineArrowColor = True
      vector         =  wrf_vector(a,wks,u_plane, v_plane,opts_vct)
  


   ;dbz
    ; optsd = res
    ;optsd@FieldTitle = "Reflectivity"
    ; optsd@cnFillPalette = "radar"
    ; optsd@cnFillOn = True
    ; optsd@cnMonoLineColor = False
    ; ;optsd@cnLineColors = "radar"
    ; res@cnLineLabelPlacementMode = "false"
    ; res@cnLineLabelDensityF=0.001
    ; optsd@lbLabelBarOn = True
    ; optsd@lbOrientation = "Vertical"
    ; optsd@pmLabelBarOrthogonalPosF= -1.25
    ; optsd@lbLabelPosition = "right"

     ;optsd@cnLineThicknessF = 1.5
     ;optsd@cnConstFEnableFill = True
     ;optsd@ContourParameters = (/ 5., 75., 5./)
     ;optsd@cnFillPalette = "radar"
     ;optsd@vcLevelPalette = "radar"
     ;optsd@cnFillColors = (/-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/)
     ;optsd@lbLabelBarOn = True
     ;optsd@gsnPanelLabelBar = True
    ; optsd@lbOrientation = "Vertical"
     ;optsd@pmLabelBarOrthogonalPosF= .05
     ;optsd@lbLabelPosition = "right"
     ;contour_dbz = wrf_contour(a,wks,dbz(1,:,:),optsd)     ; plot only lowest level

   ;ter 
    opts_ter = res
    opts_ter@FieldTitle = "Terrain Height"
    opts_ter@cnFillOn = False
    opts_ter@cnMonoLineColor = True
    ;opts_ter@cnConstFEnableFill = True
    opts_ter@cnLineThicknessF = 1.5
    opts_ter@cnLineColors = "MPL_terrain"
    opts_ter@ContourParameters = (/ 0., 250., 50./)
    contour_ter = wrf_contour(a,wks,ter,opts_ter)
    contour_ter@lbLabelBarOn = False

    plot = wrf_overlays(a,wks,(/contour_w, vector/),pltres)
end do
end
