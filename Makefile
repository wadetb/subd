#URL = file:///Users/wadeb/Documents/SubD/index.html
URL = http://wadeb.com/subd/

default:
	scp index.html root@wadeb.com:/var/www/wadeb.com/htdocs/subd/
	osascript refresh.scpt $(URL)
