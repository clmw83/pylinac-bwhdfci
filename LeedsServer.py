# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 10:19:42 2018

@author: clw34
"""

#Helpful resource: https://gist.github.com/dAnjou/2874714

import os
import matplotlib
matplotlib.use('agg')
from flask import Flask, request, redirect, url_for,make_response,render_template,Markup
from werkzeug import secure_filename
import BWHLeeds
import BWHIsoCube
import logging
import tempfile
import traceback
import base64
import io,urllib,shutil

UPLOAD_FOLDER = os.path.join('.','LeedsTemp')
if not os.path.exists(UPLOAD_FOLDER):
    os.makedirs(UPLOAD_FOLDER)
    
ALLOWED_EXTENSIONS = set(['dcm'])

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS

@app.route("/")
def index():
    return render_template("Home.html")

@app.route("/Leeds")
def leeds():
    return render_template('Leeds.html')
    
@app.route("/Leeds/Process", methods=['POST'])
def processLeeds():
    logtext=""
    if 'file' in request.files:   
        file = request.files['file']
        if file and allowed_file(file.filename):
            tempdir=tempfile.mkdtemp(dir=app.config['UPLOAD_FOLDER'])
            filename = secure_filename(file.filename)
            savename=os.path.join(tempdir, filename)
            file.save(savename)
            try:
                leeds=BWHLeeds.BWHLeeds(savename)
                leeds.analyze()
                pdfname = os.path.join(tempdir,"Leeds.pdf")
                leeds.publish_pdf(pdfname)      
                data_uri = base64.b64encode(open(pdfname, "rb").read()).decode('ascii')
                os.remove(pdfname)
                leeds.plot_analyzed_image()
                fig = matplotlib.pyplot.gcf()
                fig.set_size_inches(18.5, 10.5)
                fig.tight_layout()
                buf = io.BytesIO()
                matplotlib.pyplot.savefig(buf, format='png')
                matplotlib.pyplot.close('all')
                buf.seek(0)
                pnguri = 'data:image/png;base64,' + urllib.parse.quote(base64.b64encode(buf.read()))     
                logtext+="Successfully processed!\n"
                logtext+="Copy the below lines into the excel workbook:\n"
                logtext+=leeds.excel_lines()
                logtext+="\n"
                logtext+='<button type="button" class="btn btn-success" onclick=DownloadURI("%s")>Download PDF</button>\n'%data_uri
                logtext+='<img src = "%s" class="img-responsive"/>'%pnguri
            except:
                logtext+="Failed processing leeds!\n"
                logtext+="Are you sure you selected a Leeds phantom image?\n"
                logtext+="\nDebug information:\n"
                logtext+=traceback.format_exc()
            finally:
                os.remove(savename)                
                os.rmdir(tempdir)
        else:
            logtext="Invalid File!\n%s is not a valid DICOM file"%file.filename
    else:
        logtext="Error receiving file!\n"
    return make_response(Markup(logtext.replace('\n','<br/>')))

@app.route("/IsoCube")
def IsoCube():
    return render_template('IsoCube.html')

@app.route("/IsoCube/Process", methods=['POST'])
def processIsoCube():
    logtext=""
    uploaded_files = request.files.getlist("file")
    cubeset=BWHIsoCube.IsoCubeSet()
    
    tempdir=tempfile.mkdtemp(dir=app.config['UPLOAD_FOLDER'])
    tempfiles=[]
    if len(uploaded_files) != 8:
        logtext+="<strong>WARNING: YOU NEED TO SELECT AND UPLOAD 8 FILES (4 MV, 4KV)</strong>\n"
    logtext+="Read files:\n"
    try:
        for fs in uploaded_files:
            logtext+="%s ->"%fs.filename
            filename = secure_filename(fs.filename)
            savename=os.path.join(tempdir, filename)
            fs.save(savename)
            tempfiles.append(savename)
            etype,ang=cubeset.add_image(savename,ret=True)
            logtext+=" %s %i\n"%(etype,ang)
            os.remove(savename)
    except:
        logtext+=traceback.format_exc()
    finally:
        shutil.rmtree(tempdir)
        
    logtext+="\n"
    if len(uploaded_files) != 8:
        logtext+="<strong>The below is INVALID</strong>\n"
    logtext+="Copy into Excel:\n\n"  
    logtext+=cubeset.excel_lines()
    
    cubeset.plot()
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(10, 6)
    fig.tight_layout()
    buf = io.BytesIO()
    matplotlib.pyplot.savefig(buf, format='png')
    matplotlib.pyplot.close('all')
    buf.seek(0)
    pnguri = 'data:image/png;base64,' + urllib.parse.quote(base64.b64encode(buf.read()))     
    logtext+='<img src = "%s" class="img-responsive"/>'%pnguri

    
    return make_response(Markup(logtext.replace('\n','<br/>')))


if __name__ == "__main__":
    app.run(host="0.0.0.0",port=5001, debug=True,threaded=True)