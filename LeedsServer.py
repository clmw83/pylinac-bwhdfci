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
import logging
import tempfile
import traceback
import base64
import io,urllib

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

@app.route("/Leeds", methods=['GET', 'POST'])
def leeds():
    logtext=""
    if request.method == 'POST':   
        file = request.files['file']
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            tempdir=tempfile.mkdtemp(dir=app.config['UPLOAD_FOLDER'])
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
            return Markup(logtext.replace('\n','<br/>'))
    else:
        return render_template('Leeds.html', logtext=Markup(logtext.replace('\n','<br/>')))

@app.route("/IsoCube", methods=['GET', 'POST'])
def IsoCube():
    return render_template('IsoCube.html')

if __name__ == "__main__":
    app.run(host="0.0.0.0",port=5001, debug=True)