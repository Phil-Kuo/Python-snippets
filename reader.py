# -*- coding:utf-8 -*-

from http.cookiejar import LWPCookieJar
from urllib import request, parse
import re

#登录页面网址、验证码网址、添加名单网址
url_login       = "http://www.maxlaw.cn/Sys2013Login/chklogin.asp?Action=chk"
url_checkcode   = "http://www.maxlaw.cn/Sys2013Login/VerifyCode.asp"
url_list        = r'http://27.159.76.243:808/zMainSales/zWaiBao_Add.asp'

#cookie的相关处理和urlopener的安装
cookiejar       = LWPCookieJar()
urlopener       = request.build_opener(request.HTTPCookieProcessor(cookiejar))
request.install_opener(urlopener)
#print (cookiejar)

#模拟登录系统
request_headers = {
    'Accept-Encoding':'gzip,deflate,sdch',
    'Connection':'keep-alive',
    'Host':'www.maxlaw.cn',
    'Origin':'http://www.maxlaw.cn',
    'Referer':'http://www.maxlaw.cn/Sys2013Login/',
    'User-Agent':'Mozilla/5.0 (Windows NT 6.3; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/32.0.1700.107 Safari/537.36'
    }
post_data       = {
    'zpusername':'guoym',
    'userpass':'guoym252626',
    #'zCheckCode':checkcode,
    }

#获取验证码
response    = urlopener.open(request.Request(url_checkcode))
image       = open('check.jpg', 'wb')
image.write(response.read())
image.close()
checkcode   = None #input('Please enter the authcode:')
post_data['zCheckCode'] = checkcode
#print (cookiejar)

#发送数据
data  = parse.urlencode(post_data).encode()
res   = request.urlopen(request.Request(url_login, data, request_headers))
#print (cookiejar)

#检验登录成功
#print(res.read().decode('gbk'))
if res.getcode() == 200:
    print ("Login ok!!")
#print('***********************************************************************************')
'''
r = urlopener.open('http://120.36.188.182:808/zMainSales/zWaiBao_Add.asp')
#print (cookiejar)
print(r.geturl())
print(r.read(1000).decode('gbk'))
#print('***********************************************************************************')
'''
def list_add(tel_no, name):
    #global url_list
    list_headers = {
        'Host':'27.159.76.243:808',
        'Origin':'http://27.159.76.243:808',
        'Referer':'http://27.159.76.243:808/zMainSales/zWaiBao_Add.asp',
        'Cookie':'ASPSESSIONIDACSSRDBC=PHMLNOOAAKCIGGLIEONNKBHH',  ##此处的cookie是可以添加表单的关键，并且可能每天都不太一样，但可能可以持续一天不变。
        'User-Agent':'Mozilla/5.0 (Windows NT 6.3; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/36.0.1985.143 Safari/537.36',
        }

    list_data    = {
        'zSaveTel':tel_no,
        'zSaveName':name,
        'zSalesID':'685451259973',
        'zSalesKey':'guoym',
        'action1':'save',
        }
    datas  = parse.urlencode(list_data).encode()
    response   = request.urlopen(request.Request(url_list, datas, headers = list_headers))
    #print (cookiejar)
    message = response.read(20).decode('gbk')
    #print(message)
    return message
    
#添加名单
import csv
with open ('new.csv', 'r') as f:
    reader = csv.reader(f)
    i = 0
    add_list_failed = 0
    for row in reader:
        #print(row[0])
        tel_no = row[1]
        name_1   = row[0]
        name   = name_1.encode('gbk')
        print(tel_no,name_1)
        msg = list_add(tel_no, name)
        try:
            if re.search('Error', msg).group() == 'Error':
                add_list_failed += 1
                print('error')
        except AttributeError:
            pass
        i += 1
        
print ("total:", i, "failed:", add_list_failed)






